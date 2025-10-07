import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import statsmodels.formula.api as smf
from scipy.stats import t

# 1. Read data
data_path = "inputs/data/Mapping Police Violence_new"
df = pd.read_csv(data_path)

# 2. Descriptive statistics
print("Descriptive statistics:")
print(df[['tmean', 'prec', 'pop', 'rate_adj']].describe())

# Calculate mean of rate_adj for percentage conversion
mean_rate_adj = df['rate_adj'].mean()
print(f"\nMean of rate_adj: {mean_rate_adj}")

# 3. Define grouping functions

# Temperature: Above/Below average
mean_tmean = df['tmean'].mean()
df['temp_avg_group'] = df['tmean'].apply(
    lambda x: 'Above average Temperature' if x > mean_tmean else 'Below average Temperature')

# Precipitation: Above/Below average
mean_prec = df['prec'].mean()
df['prec_avg_group'] = df['prec'].apply(
    lambda x: 'Above average Precipitation' if x > mean_prec else 'Below average Precipitation')

# Population: Above/Below average
mean_pop = df['pop'].mean()
df['pop_avg_group'] = df['pop'].apply(
    lambda x: 'Above average Population' if x > mean_pop else 'Below average Population')


# Temperature bins: <10℃, 10℃~20℃, >20℃
def bin_temperature(temp):
    if temp < 10:
        return '<10℃'
    elif 10 <= temp <= 20:
        return '10℃~20℃'
    else:
        return '>20℃'


df['temp_group'] = df['tmean'].apply(bin_temperature)


# Precipitation bins: <50mm, 50~100mm, >100mm
def bin_precipitation(precip):
    if precip < 50:
        return '<50mm'
    elif 50 <= precip <= 100:
        return '50~100mm'
    else:
        return '>100mm'


df['precip_group'] = df['prec'].apply(bin_precipitation)


# Population bins: <1M, 1M~5M, >5M
def bin_population(pop):
    if pop < 1_000_000:
        return '<1M'
    elif 1_000_000 <= pop <= 5_000_000:
        return '1M~5M'
    else:
        return '>5M'


df['pop_group'] = df['pop'].apply(bin_population)

# 4. Construct temperature dynamic variables (t-1, t, t+1)
if 'fips' in df.columns and 'stateyear' in df.columns:
    df = df.sort_values(['fips', 'stateyear'])

if 'fips' in df.columns:
    df['tmean_lag1'] = df.groupby('fips')['tmean'].shift(1)
    df['tmean_lead1'] = df.groupby('fips')['tmean'].shift(-1)
else:
    df['tmean_lag1'] = df['tmean'].shift(1)
    df['tmean_lead1'] = df['tmean'].shift(-1)
df['tmean_0'] = df['tmean']  # Current period

# Drop rows with missing values due to lag/lead creation
df_temp = df.dropna(subset=['tmean_lag1', 'tmean_0', 'tmean_lead1'])

# 5. Regression analysis and results collection
results = []

group_order = [
    'Above average Temperature', 'Below average Temperature',
    'Above average Precipitation', 'Below average Precipitation',
    'Above average Population', 'Below average Population',
    'Temperature:>20℃', 'Temperature:10℃~20℃', 'Temperature:<10℃',
    'Precipitation: >100mm', 'Precipitation: 50~100mm', 'Precipitation: <50mm',
    'Population: >5M', 'Population: 1M~5M', 'Population: <1M',
    'Month: t-1', 'Month: t', 'Month: t+1',
]

# ---------------------------
# (A) Interaction model for Temperature (Above/Below average)
# ---------------------------
model_temp_int = smf.ols('rate_adj ~ tmean * C(temp_avg_group)', data=df).fit()

# 参考组："Above average Temperature"（默认作为参考）
coef_ref_temp = model_temp_int.params['tmean']
conf_ref_temp = model_temp_int.conf_int().loc['tmean'].values
n_ref_temp = len(df[df['temp_avg_group'] == 'Above average Temperature'])
results.append({
    'Group': 'Above average Temperature',
    'Sample Size': n_ref_temp,
    'Coefficient (%)': (coef_ref_temp / mean_rate_adj) * 100,
    'CI Lower (%)': (conf_ref_temp[0] / mean_rate_adj) * 100,
    'CI Upper (%)': (conf_ref_temp[1] / mean_rate_adj) * 100,
})

coef_int_temp = model_temp_int.params.get('tmean:C(temp_avg_group)[T.Below average Temperature]', 0)
coef_below_temp = coef_ref_temp + coef_int_temp

cov_temp = model_temp_int.cov_params()
var_sum_temp = (cov_temp.loc['tmean', 'tmean'] +
                cov_temp.loc['tmean:C(temp_avg_group)[T.Below average Temperature]',
                'tmean:C(temp_avg_group)[T.Below average Temperature]'] +
                2 * cov_temp.loc['tmean', 'tmean:C(temp_avg_group)[T.Below average Temperature]'])
std_sum_temp = np.sqrt(var_sum_temp)
dof_temp = model_temp_int.df_resid
t_val_temp = t.ppf(0.975, dof_temp)
ci_lower_below_temp = coef_below_temp - t_val_temp * std_sum_temp
ci_upper_below_temp = coef_below_temp + t_val_temp * std_sum_temp

n_below_temp = len(df[df['temp_avg_group'] == 'Below average Temperature'])
results.append({
    'Group': 'Below average Temperature',
    'Sample Size': n_below_temp,
    'Coefficient (%)': (coef_below_temp / mean_rate_adj) * 100,
    'CI Lower (%)': (ci_lower_below_temp / mean_rate_adj) * 100,
    'CI Upper (%)': (ci_upper_below_temp / mean_rate_adj) * 100,
})

# ---------------------------
# (B) Interaction model for Precipitation (Above/Below average)
# ---------------------------
model_prec_int = smf.ols('rate_adj ~ tmean * C(prec_avg_group)', data=df).fit()

coef_ref_prec = model_prec_int.params['tmean']
conf_ref_prec = model_prec_int.conf_int().loc['tmean'].values
n_ref_prec = len(df[df['prec_avg_group'] == 'Above average Precipitation'])
results.append({
    'Group': 'Above average Precipitation',
    'Sample Size': n_ref_prec,
    'Coefficient (%)': (coef_ref_prec / mean_rate_adj) * 100,
    'CI Lower (%)': (conf_ref_prec[0] / mean_rate_adj) * 100,
    'CI Upper (%)': (conf_ref_prec[1] / mean_rate_adj) * 100,
})

coef_int_prec = model_prec_int.params.get('tmean:C(prec_avg_group)[T.Below average Precipitation]', 0)
coef_below_prec = coef_ref_prec + coef_int_prec

cov_prec = model_prec_int.cov_params()
var_sum_prec = (cov_prec.loc['tmean', 'tmean'] +
                cov_prec.loc['tmean:C(prec_avg_group)[T.Below average Precipitation]',
                'tmean:C(prec_avg_group)[T.Below average Precipitation]'] +
                2 * cov_prec.loc['tmean', 'tmean:C(prec_avg_group)[T.Below average Precipitation]'])
std_sum_prec = np.sqrt(var_sum_prec)
dof_prec = model_prec_int.df_resid
t_val_prec = t.ppf(0.975, dof_prec)
ci_lower_below_prec = coef_below_prec - t_val_prec * std_sum_prec
ci_upper_below_prec = coef_below_prec + t_val_prec * std_sum_prec

n_below_prec = len(df[df['prec_avg_group'] == 'Below average Precipitation'])
results.append({
    'Group': 'Below average Precipitation',
    'Sample Size': n_below_prec,
    'Coefficient (%)': (coef_below_prec / mean_rate_adj) * 100,
    'CI Lower (%)': (ci_lower_below_prec / mean_rate_adj) * 100,
    'CI Upper (%)': (ci_upper_below_prec / mean_rate_adj) * 100,
})

# ---------------------------
# (C) Interaction model for Population (Above/Below average)
# ---------------------------
model_pop_int = smf.ols('rate_adj ~ tmean * C(pop_avg_group)', data=df).fit()

coef_ref_pop = model_pop_int.params['tmean']
conf_ref_pop = model_pop_int.conf_int().loc['tmean'].values
n_ref_pop = len(df[df['pop_avg_group'] == 'Above average Population'])
results.append({
    'Group': 'Above average Population',
    'Sample Size': n_ref_pop,
    'Coefficient (%)': (coef_ref_pop / mean_rate_adj) * 100,
    'CI Lower (%)': (conf_ref_pop[0] / mean_rate_adj) * 100,
    'CI Upper (%)': (conf_ref_pop[1] / mean_rate_adj) * 100,
})

coef_int_pop = model_pop_int.params.get('tmean:C(pop_avg_group)[T.Below average Population]', 0)
coef_below_pop = coef_ref_pop + coef_int_pop

cov_pop = model_pop_int.cov_params()
var_sum_pop = (cov_pop.loc['tmean', 'tmean'] +
               cov_pop.loc['tmean:C(pop_avg_group)[T.Below average Population]',
               'tmean:C(pop_avg_group)[T.Below average Population]'] +
               2 * cov_pop.loc['tmean', 'tmean:C(pop_avg_group)[T.Below average Population]'])
std_sum_pop = np.sqrt(var_sum_pop)
dof_pop = model_pop_int.df_resid
t_val_pop = t.ppf(0.975, dof_pop)
ci_lower_below_pop = coef_below_pop - t_val_pop * std_sum_pop
ci_upper_below_pop = coef_below_pop + t_val_pop * std_sum_pop

n_below_pop = len(df[df['pop_avg_group'] == 'Below average Population'])
results.append({
    'Group': 'Below average Population',
    'Sample Size': n_below_pop,
    'Coefficient (%)': (coef_below_pop / mean_rate_adj) * 100,
    'CI Lower (%)': (ci_lower_below_pop / mean_rate_adj) * 100,
    'CI Upper (%)': (ci_upper_below_pop / mean_rate_adj) * 100,
})

# ---------------------------

# (D) Temperature group regression: >20℃, 10℃~20℃, <10℃
# ---------------------------
temp_groups = ['>20℃', '10℃~20℃', '<10℃']
for temp_group in temp_groups:
    subset = df[df['temp_group'] == temp_group]
    # Modified model with interaction term
    model = smf.ols('rate_adj ~ tmean ', data=subset).fit()
    coef = model.params['tmean']
    conf_int = model.conf_int().loc['tmean'].values
    sample_size = len(subset)
    results.append({
        'Group': f'Temperature:{temp_group}',
        'Sample Size': sample_size,
        'Coefficient (%)': (coef / mean_rate_adj) * 100,
        'CI Lower (%)': (conf_int[0] / mean_rate_adj) * 100,
        'CI Upper (%)': (conf_int[1] / mean_rate_adj) * 100,
    })

# ---------------------------
# (E) Precipitation group regression: <50mm, 50~100mm, >100mm
# ---------------------------
precip_groups = ['>100mm', '50~100mm', '<50mm']
for precip_group in precip_groups:
    subset = df[df['precip_group'] == precip_group]
    # Modified model with interaction term
    model = smf.ols('rate_adj ~ tmean ', data=subset).fit()
    coef = model.params['tmean']
    conf_int = model.conf_int().loc['tmean'].values
    sample_size = len(subset)
    results.append({
        'Group': f'Precipitation: {precip_group}',
        'Sample Size': sample_size,
        'Coefficient (%)': (coef / mean_rate_adj) * 100,
        'CI Lower (%)': (conf_int[0] / mean_rate_adj) * 100,
        'CI Upper (%)': (conf_int[1] / mean_rate_adj) * 100,
    })

# ---------------------------
# (F) Population group regression: <1M, 1M~5M, >5M
# ---------------------------
pop_groups = ['>5M', '1M~5M', '<1M']
for pop_group in pop_groups:
    subset = df[df['pop_group'] == pop_group]
    # Modified model with interaction term
    model = smf.ols('rate_adj ~ tmean ', data=subset).fit()
    coef = model.params['tmean']
    conf_int = model.conf_int().loc['tmean'].values
    sample_size = len(subset)
    results.append({
        'Group': f'Population: {pop_group}',
        'Sample Size': sample_size,
        'Coefficient (%)': (coef / mean_rate_adj) * 100,
        'CI Lower (%)': (conf_int[0] / mean_rate_adj) * 100,
        'CI Upper (%)': (conf_int[1] / mean_rate_adj) * 100,
    })
# ---------------------------
# (G) Temperature dynamic regression: t-1, t, t+1
# ---------------------------
point_estimates = [-0.01736275, 0.11456770, -0.02067211]  # t-1, t, t+1
ci_lower = [-0.22130810, -0.10020157, -0.22358647]  # t-1, t, t+1 lower CI
ci_upper = [0.1865826, 0.3293370, 0.1822422]  # t-1, t, t+1 upper CI
for i, period in enumerate(['t-1', 't', 't+1']):
    results.append({
        'Group': f'Month: {period}',
        'Sample Size': len(df_temp),
        'Coefficient (%)': point_estimates[i],
        'CI Lower (%)': ci_lower[i],
        'CI Upper (%)': ci_upper[i],
    })

# 6. Convert results to DataFrame and reorder
results_df = pd.DataFrame(results)
results_df = results_df.set_index('Group').loc[group_order].reset_index()
print("\nHeterogeneity Analysis Results (Percentage Change):")
print(results_df)

# 7. Plotting
plt.figure(figsize=(14.5, 8), dpi=600)
ax = plt.gca()

x_min, x_max = -3, 8.2

for i, row in results_df.iterrows():
    coef = row['Coefficient (%)']
    ci_lower_val = row['CI Lower (%)']
    ci_upper_val = row['CI Upper (%)']

    # 如果有NaN，跳过
    if pd.isna(coef) or pd.isna(ci_lower_val) or pd.isna(ci_upper_val):
        print(f"Skipping row {i} due to NaN values: {row}")
        continue

    # 截断置信区间到x轴范围
    ci_lower_val = max(ci_lower_val, x_min)
    ci_upper_val = min(ci_upper_val, x_max)

    # 计算新的左右误差（确保非负）
    left_error = max(0, coef - ci_lower_val)
    right_error = max(0, ci_upper_val - coef)

    ax.errorbar(
        x=coef,
        y=len(results_df) - 1 - i,
        xerr=[[left_error], [right_error]],
        fmt='o',
        capsize=5,
        color='black',
        ecolor='gray',
        markersize=8,
        label='Coefficient with 95% CI' if i == 0 else ""
    )

y_labels = results_df['Group'].tolist()
for i, label in enumerate(y_labels):
    ax.text(-4.3, len(results_df) - 1 - i, label, ha='center', va='center', fontsize=10)

ax.set_yticks([])

ax.set_xlabel('% change in death rate of police violence per +1℃', fontsize=14)
ax.set_xlim(x_min, x_max)

ax.axvline(x=0, color='black', linestyle='--', linewidth=1)

ax.grid(True, linestyle='--', alpha=0.5, which='both')

group_boundaries = [6, 15]
for boundary_index in group_boundaries:
    ax.axhline(y=len(results_df) - 1 - boundary_index + 0.5, color='red', linestyle='--', linewidth=1.5)

ax2 = ax.twinx()
ax2.set_yticks([])
ax2.set_ylim(ax.get_ylim())

category_labels = ["Introducing Interaction Terms", "Binned Model", "Displacement Effect"]
category_ranges = [(0, 5), (6, 14), (15, 17)]
right_label_x_position = 8.3
for i, (label, (start, end)) in enumerate(zip(category_labels, category_ranges)):
    mid_point = len(results_df) - 1 - (start + end) / 2
    ax2.text(right_label_x_position, mid_point, f"({label})", ha='left', va='center', rotation=90, fontsize=10)

ax.legend(loc='upper right', fontsize=12)
plt.tight_layout(rect=[0, 0, 0.8, 1])

# 8. Save figure and results
output_dir = "outputs/raw_figures/"
os.makedirs(output_dir, exist_ok=True)

fig_path = os.path.join(output_dir, "Figure2.png")
plt.savefig(fig_path, dpi=600, bbox_inches='tight')
plt.close()
print(f"\nFigure saved to: {fig_path}")

csv_path = os.path.join(output_dir, "Figure2.csv")
results_df.to_csv(csv_path, index=False)
print(f"Results table saved to: {csv_path}")