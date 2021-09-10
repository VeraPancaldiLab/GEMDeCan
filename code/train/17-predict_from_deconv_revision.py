# %%
"""
Train ElasticNet penalized Logistic regression on deconvolution data
only to assess how important/ usefull new signatures are.
"""

# %%
import pandas as pd
import numpy as np

import random
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from pathlib import Path
from time import time
from multiprocessing import cpu_count
import os
import joblib
from time import time

import xgboost as xgb
from sklearn import tree
from sklearn import ensemble
from sklearn import svm
# from sklearn.linear_model import lasso_path, enet_path
from sklearn import linear_model

from sklearn.preprocessing import StandardScaler, PowerTransformer, RobustScaler
from sklearn.model_selection import train_test_split, cross_val_score
import sklearn.metrics as metrics

from sklearn.experimental import enable_halving_search_cv  # noqa
from sklearn.model_selection import cross_validate, GridSearchCV, RandomizedSearchCV, HalvingGridSearchCV
from sklearn.utils.class_weight import compute_sample_weight
from scipy.stats import loguniform

from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer, SimpleImputer, KNNImputer
# %%
MERGE_OLD_SIG = False
ONLY_MELANOMA = True
ALL_SIG = True
DROP_SCORES = True   # discard XCELL's 'ImmuneScore', 'StromaScore' and 'MicroenvironmentScore'

path_interm = Path("../../data/intermediate/Deconvolution_paper")
if MERGE_OLD_SIG:
    dir_save = Path("../../data/processed/Deconvolution_paper_revision/new_and_old_signatures")
elif ALL_SIG:
    dir_save = Path("../../data/processed/Deconvolution_paper_revision/all_signatures")
if ONLY_MELANOMA:
    dir_save = dir_save / "only_melanoma"
else:
    dir_save = dir_save / "all_cancer_types"
if DROP_SCORES:
    dir_save = dir_save / "drop_scores"
else:
    dir_save = dir_save / "keep_scores"
    
dir_save.mkdir(parents=True, exist_ok=True)

response = pd.read_csv(path_interm / "response_train.csv", index_col=0)

datasets = ['Gide', 'Riaz', 'Hugo', 'Cloughesy']
deconv = None
for dataset in datasets:
    path_data = path_interm / 'revision_deconv' / ('all_deconvolutions_' + dataset + '.txt')
    if deconv is None:
        deconv = pd.read_csv(path_data, sep='\t', index_col=0)
    else:
        new_deconv = pd.read_csv(path_data, sep='\t', index_col=0)
        deconv = deconv.append(new_deconv)



common_ids = set(response.index).intersection(deconv.index)
df_all = deconv.loc[common_ids, :]
response = response.loc[df_all.index, :]

# drop variables with unique values
drop_col = [x for x in df_all.columns if df_all[x].unique().size == 1]
print("dropping these columns because they have a unique value:")
for i in drop_col:
    print("    ", i)
matches = [] #['XCELL', 'TAMs', 'CD226']
if DROP_SCORES:
    matches = matches + ['Score']
drop_col = drop_col + [x for x in df_all.columns if any(y in x for y in matches)]
if len(drop_col) > 0:
    df_all.drop(columns=drop_col, inplace=True)
# drop Maha samples
drop_row = [x for x in df_all.index if 'Maha' in x]
if ONLY_MELANOMA:
    drop_row = drop_row + [x for x in df_all.index if 'Snyder' in x]
    drop_row = drop_row + [x for x in df_all.index if 'Cloughesy' in x]
df_all.drop(labels=drop_row, axis=0, inplace=True)
response.drop(labels=drop_row, inplace=True)
# set inf values to nan
to_nan = ~np.isfinite(df_all.values)
nb_nan = to_nan.sum()
if nb_nan != 0:
    print(f"There are {nb_nan} nan values")
    df_all[to_nan] = np.nan
    # impute non finite values (nan, +/-inf)
    # imputer = IterativeImputer(max_iter=100, random_state=0)
    imputer = KNNImputer(n_neighbors=5, weights="distance")
    df_all.loc[:,:] = imputer.fit_transform(df_all.values)
common_col = df_all.columns.values
nb_var = df_all.shape[1]

# rename variables
dic_rename = {
    'Epidish': 'EpiDISH',
    'DeconRNASeq': 'deconRNASeq',
    'Quantiseq': 'quanTIseq',
}
for key, val in dic_rename.items():
    new_cols = [x.replace(key, val) for x in df_all.columns]
    df_all.columns = new_cols
# %%
y_pred

# %%
condi = conditions[0]
test_name = 'Gide'

var_idx = [x for x in df_all.columns if x.startswith(condi)]
str_condi = condi.strip('_')

score_split = {}


print(f"    LODO training using {test_name} as test set")
test_ids = [x for x in df_all.index if test_name in x]
train_ids = [x for x in df_all.index if test_name not in x]
X_train = df_all.loc[train_ids, var_idx].values
X_test= df_all.loc[test_ids, var_idx].values
y_train = response['BOR - binary'].loc[train_ids].values
y_test = response['BOR - binary'].loc[test_ids].values
# reverse labels to match "equaling progressive disease" objective
y_train = -(y_train-1)
y_test = -(y_test-1)
# Standardize data to give same weight to regularization
# scaler = PowerTransformer(method='yeo-johnson')
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

clf = linear_model.LogisticRegressionCV(
    Cs=20, 
    penalty='elasticnet', 
    # scoring='neg_log_loss', 
    scoring='roc_auc', 
    solver='saga', 
    l1_ratios=l1_ratios,
    max_iter=10000,
    n_jobs=-1,  # or n_jobs-1 to leave one core available
)
clf = clf.fit(X_train, y_train)

y_pred = clf.predict_proba(X_test)[:, 1]
score = metrics.roc_auc_score(y_test, y_pred)
score_split[test_name] = score
print(f"         ROC AUC {test_name}: {score}")

# Save model coefficients and plots
param_string = f"signature-{str_condi}_split-lodo-{test_name}"
l1_ratio = np.round(clf.l1_ratio_[0], decimals=4)
C = np.round(clf.C_[0], decimals=4)

coef = pd.DataFrame({'coef': clf.coef_.flatten()}, index=var_idx)
coef['abs coef'] = coef['coef'].abs()
coef = coef.sort_values(by='abs coef', ascending=False)
coef['% total'] = coef['abs coef'] / coef['abs coef'].sum()
coef['cum % total'] = coef['% total'].cumsum()
nb_coef = coef.shape[0]

nb_coef_plot = min(20, nb_coef)
labels = coef.index[:nb_coef_plot]
plt.figure()
ax = coef.loc[labels, 'cum % total'].plot()
ticks_pos = np.linspace(start=0, stop=nb_coef_plot-1, num=nb_coef_plot)
# ticks_label = np.round(ticks_label, decimals=2)
ax.set_xticks(ticks_pos)
# ax.set_xticklabels(ticks_label)
ax.set_xticklabels(labels, rotation=45, ha='right')
plt.xlabel('variables')
plt.ylabel('cumulative % coef')
plt.title(param_string + f" l1_ratio {l1_ratio}, C {C}")

# %% Compare models trained on single signatures

if ALL_SIG:
    # we use trailing underscores to avoid including derived signatures
    # like EpiDISH_BPRNACan3Dprom --> EpiDISH_BPRNACan3Dprom-enhan
    conditions = [
        'quanTIseq',
        'MCP',
        'XCELL',
        'EpiDISH_BPRNACan_',
        'deconRNASeq_BPRNACan_',
        'EpiDISH_BPRNACanProMet_',
        'deconRNASeq_BPRNACanProMet_',
        'EpiDISH_BPRNACan3DProMet_',
        'deconRNASeq_BPRNACan3DProMet_',
        'EpiDISH_CBSX_CIBERSORTx_ref_rna-seq_',
        'deconRNASeq_CBSX_CIBERSORTx_ref_rna-seq_',
        'EpiDISH_CBSX_NSCLC_single-cell_',
        'deconRNASeq_CBSX_NSCLC_single-cell_',
        'CBSX_CIBERSORTx_ref_rna-seq_',
        'CBSX_NSCLC_single-cell_',
    ]
else:
    conditions = [
        'quanTIseq',
        'MCP',
        'EpiDISH_BPprom_314',
        'DeconRNA_BPprom_314',
        'EpiDISH_BPRNACan',
        'DeconRNA_BPRNACan',
        'EpiDISH_BPRNACan3DMet',
        'DeconRNA_BPRNACan3DMet',
        'EpiDISH_BPRNACanProMet3D',
        'DeconRNA_BPRNACanProMet3D',
    ]

if ONLY_MELANOMA or MERGE_OLD_SIG:
    datasets = ['Gide', 'Hugo', 'Riaz']
else:
    datasets = ['Gide', 'Hugo', 'Riaz', 'Snyder']
# l1_ratio = 0 the penalty is an L2 penalty [Ridge]
# l1_ratio = 1 the penalty is an L1 penalty (Lasso)

# Test either one of those combinations
# default parameter
l1_ratios_list = [
    ['default', [0.5]],
    ['naive', np.linspace(0, 1, 21)],           # naive param grid
    ['advised', [.1, .5, .7, .9, .95, .99, 1]], # advised in scikit-learn documentation
]

start = time()

dir_save_orig = dir_save
# l1_name, l1_ratios = l1_ratios_list[0]
for l1_name, l1_ratios in l1_ratios_list:
    dir_save = dir_save_orig / f'l1_ratios-{l1_name}'
    dir_save.mkdir(parents=True, exist_ok=True)

    score_condi = {}
    # condi = conditions[0]
    for condi in conditions:
        print(f"Training using {condi}")
        var_idx = [x for x in df_all.columns if x.startswith(condi)]
        str_condi = condi.strip('_')
        
        score_split = {}
        # LODO splitting
        # test_name = datasets[0]
        for test_name in datasets:
            print(f"    LODO training using {test_name} as test set")
            test_ids = [x for x in df_all.index if test_name in x]
            train_ids = [x for x in df_all.index if test_name not in x]
            X_train = df_all.loc[train_ids, var_idx].values
            X_test= df_all.loc[test_ids, var_idx].values
            y_train = response['BOR - binary'].loc[train_ids].values
            y_test = response['BOR - binary'].loc[test_ids].values
            # reverse labels to match "equaling progressive disease" objective
            y_train = -(y_train-1)
            y_test = -(y_test-1)
            # Standardize data to give same weight to regularization
            # scaler = PowerTransformer(method='yeo-johnson')
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)

            clf = linear_model.LogisticRegressionCV(
                Cs=20, 
                penalty='elasticnet', 
                # scoring='neg_log_loss', 
                scoring='roc_auc', 
                solver='saga', 
                l1_ratios=l1_ratios,
                max_iter=10000,
                n_jobs=-1,  # or n_jobs-1 to leave one core available
            )
            clf = clf.fit(X_train, y_train)

            y_pred = clf.predict_proba(X_test)[:, 1]
            score = metrics.roc_auc_score(y_test, y_pred)
            score_split[test_name] = score
            print(f"         ROC AUC {test_name}: {score}")

            # Save model coefficients and plots
            param_string = f"signature-{str_condi}_split-lodo-{test_name}"
            l1_ratio = np.round(clf.l1_ratio_[0], decimals=4)
            C = np.round(clf.C_[0], decimals=4)

            coef = pd.DataFrame({'coef': clf.coef_.flatten()}, index=var_idx)
            coef['abs coef'] = coef['coef'].abs()
            coef = coef.sort_values(by='abs coef', ascending=False)
            coef['% total'] = coef['abs coef'] / coef['abs coef'].sum()
            coef['cum % total'] = coef['% total'].cumsum()
            coef.to_csv(dir_save / f"LogisticRegressionCV_coefficients_{param_string}.csv")
            nb_coef = coef.shape[0]

            nb_coef_plot = min(20, nb_coef)
            labels = coef.index[:nb_coef_plot]
            plt.figure()
            ax = coef.loc[labels, 'cum % total'].plot()
            ticks_pos = np.linspace(start=0, stop=nb_coef_plot-1, num=nb_coef_plot)
            # ticks_label = np.round(ticks_label, decimals=2)
            ax.set_xticks(ticks_pos)
            # ax.set_xticklabels(ticks_label)
            ax.set_xticklabels(labels, rotation=45, ha='right')
            plt.xlabel('variables')
            plt.ylabel('cumulative % coef')
            plt.title(param_string + f" l1_ratio {l1_ratio}, C {C}")
            plt.savefig(dir_save / f"cumulative_abs_coef_{param_string}.png", bbox_inches='tight', facecolor='white')

            plt.figure()
            ax = coef.loc[labels, 'coef'].plot()
            ax.hlines(y=0, xmin=0, xmax=nb_coef_plot-1, colors='gray', linestyles='dashed')
            ticks_pos = np.linspace(start=0, stop=nb_coef_plot-1, num=nb_coef_plot)
            # ticks_label = np.round(ticks_label, decimals=2)
            ax.set_xticks(ticks_pos)
            # ax.set_xticklabels(ticks_label)
            ax.set_xticklabels(labels, rotation=45, ha='right')
            plt.xlabel('variables')
            plt.ylabel('coef')
            plt.title(param_string + f" l1_ratio {l1_ratio}, C {C}")
            plt.savefig(dir_save / f"coef_{param_string}.png", bbox_inches='tight', facecolor='white')


        # ------ All datasets with CV splitting ------
        print(f"    CV training on all datasets")
        X = df_all.loc[:, var_idx].values
        y = response['BOR - binary'].values
        # stratify train / test by dataset and response
        strat = response.apply(lambda x:x['dataset'] + '-' + str(x['BOR - binary']), axis=1)
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, 
            test_size=0.25, 
            random_state=0, 
            shuffle=True, 
            stratify=strat,
        )
        # reverse labels to match "equaling progressive disease" objective
        y_train = -(y_train-1)
        y_test = -(y_test-1)
        # Standardize data to give same weight to regularization
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        clf = linear_model.LogisticRegressionCV(
            Cs=20, 
            penalty='elasticnet', 
            # scoring='neg_log_loss', 
            scoring='roc_auc', 
            solver='saga', 
            l1_ratios=l1_ratios,
            max_iter=10000,
            n_jobs=-1,  # or n_jobs-1 to leave one core available
        )
        clf = clf.fit(X_train, y_train)

        y_pred = clf.predict_proba(X_test)[:, 1]
        score = metrics.roc_auc_score(y_test, y_pred)
        score_split['CV_all_datasets'] = score
        print(f"         ROC AUC CV all datasets: {score}")
        score_condi[str_condi] = score_split

        # Save model coefficients and plots
        param_string = f"signature-{str_condi}_split-CV-5folds"
        l1_ratio = np.round(clf.l1_ratio_[0], decimals=4)
        C = np.round(clf.C_[0], decimals=4)

        coef = pd.DataFrame({'coef': clf.coef_.flatten()}, index=var_idx)
        coef['abs coef'] = coef['coef'].abs()
        coef = coef.sort_values(by='abs coef', ascending=False)
        coef['% total'] = coef['abs coef'] / coef['abs coef'].sum()
        coef['cum % total'] = coef['% total'].cumsum()
        coef.to_csv(dir_save / f"LogisticRegressionCV_coefficients_{param_string}.csv")
        nb_coef = coef.shape[0]

        nb_coef_plot = min(20, nb_coef)
        labels = coef.index[:nb_coef_plot]
        plt.figure()
        ax = coef.loc[labels, 'cum % total'].plot()
        ticks_pos = np.linspace(start=0, stop=nb_coef_plot-1, num=nb_coef_plot)
        # ticks_label = np.round(ticks_label, decimals=2)
        ax.set_xticks(ticks_pos)
        # ax.set_xticklabels(ticks_label)
        ax.set_xticklabels(labels, rotation=45, ha='right')
        plt.xlabel('variables')
        plt.ylabel('cumulative % coef')
        plt.title(param_string + f" l1_ratio {l1_ratio}, C {C}")
        plt.savefig(dir_save / f"cumulative_abs_coef_{param_string}.png", bbox_inches='tight', facecolor='white')

        plt.figure()
        ax = coef.loc[labels, 'coef'].plot()
        ax.hlines(y=0, xmin=0, xmax=nb_coef_plot-1, colors='gray', linestyles='dashed')
        ticks_pos = np.linspace(start=0, stop=nb_coef_plot-1, num=nb_coef_plot)
        # ticks_label = np.round(ticks_label, decimals=2)
        ax.set_xticks(ticks_pos)
        # ax.set_xticklabels(ticks_label)
        ax.set_xticklabels(labels, rotation=45, ha='right')
        plt.xlabel('variables')
        plt.ylabel('coef')
        plt.title(param_string + f" l1_ratio {l1_ratio}, C {C}")
        plt.savefig(dir_save / f"coef_{param_string}.png", bbox_inches='tight', facecolor='white')

        plt.show()
        plt.close('all')

    scores = pd.DataFrame.from_dict(score_condi, orient='index')
    scores.index.name = 'signature'
    scores.to_csv(dir_save / 'ROC_AUC_single-signature.csv')

end = time()
duration = end - start
print("\n"*5)
print("-------------------------------------------")
print("------------------- END -------------------")
print("-------------------------------------------")
print(f"Training took {duration}s")
# %% Check cell types proportions
# from https://stackoverflow.com/questions/45664519/export-pandas-styled-table-to-image-file
# install https://wkhtmltopdf.org/ 
# and  https://github.com/jarrekk/imgkit

import imgkit

dir_save = dir_save_orig / f'cell_proportions'
dir_save.mkdir(parents=True, exist_ok=True)
for condi in conditions:
    var_idx = [x for x in df_all.columns if condi in x]
    str_condi = condi.strip('_')
    prop = df_all.loc[:, var_idx].mean(axis=0).to_frame(name='all')
    # LODO splitting
    # test_name = datasets[0]
    for test_name in datasets:
        test_ids = [x for x in df_all.index if test_name in x]
        prop = prop.join(df_all.loc[test_ids, var_idx].mean(axis=0).to_frame(name=test_name))
    prop.to_csv(dir_save / f'proportions_signature-{str_condi}.csv')
    styled_table = prop.style.background_gradient(cmap='RdYlGn_r', axis=1)
    html = styled_table.render()
    imgkit.from_string(html, dir_save / f'proportions_signature-{str_condi}.png')

# %% Targeted comparison
# condi1 = 'EpiDISH_BPRNACan3Dprom'
# condi2 = 'deconRNASeq_BPRNACan'
condi1 = 'EpiDISH_BPRNACan3DProMet'
condi2 = 'deconRNASeq_BPRNACan'
prop1 = pd.read_csv(dir_save / f'proportions_signature-{condi1}.csv', index_col=0)[['all']]
prop2 = pd.read_csv(dir_save / f'proportions_signature-{condi2}.csv', index_col=0)[['all']]

new_index = [x.split('_')[-1] for x in prop1.index]
prop1.index = new_index
prop1.columns = [condi1]
new_index = [x.split('_')[-1] for x in prop2.index]
prop2.index = new_index
prop2.columns = [condi2]

prop = prop1.join(prop2)
str_condi = '-'.join([condi1, condi2])
prop.to_csv(dir_save / f'proportions_signature-{str_condi}.csv')
styled_table = prop.style.background_gradient(cmap='RdYlGn_r', axis=1)
html = styled_table.render()
imgkit.from_string(html, dir_save / f'proportions_signature-{str_condi}.png')

# %%

coef_1 = pd.read_csv("../../data/processed/Deconvolution_paper_revision/all_signatures/no_bladder/l1_ratios-default/LogisticRegressionCV_coefficients_signature-deconRNASeq_BPRNACan3DProMet_split-CV-5folds.csv", index_col=0)
new_index = [x.replace('deconRNASeq_BPRNACan3DProMet_', '') for x in coef_1.index]
# new_index = [x.replace('BPRNACan3Dprom', 'BPRNACanProMet') for x in coef.index]
# new_index = [x.replace('DeconRNA', 'deconRNAseq') for x in coef.index]
coef_1.index = new_index
nb_coef = coef_1.shape[0]
coef_1.index.name = 'cell type'

coef_2 = pd.read_csv("../../data/processed/Deconvolution_paper_revision/all_signatures/no_bladder/l1_ratios-default/LogisticRegressionCV_coefficients_signature-deconRNAseq_BPRNACan_split-CV-5folds.csv", index_col=0)
new_index = [x.replace('deconRNAseq_BPRNACan_', '') for x in coef_2.index]
coef_2.index = new_index
coef_2.index.name = 'cell type'

# make same order of variable for the 2 plots, by max abs value
abs_coef = coef_1['abs coef'] + coef_2['abs coef']
abs_coef.sort_values(ascending=False, inplace=True)
coef_1 = coef_1.loc[abs_coef.index, :]
coef_2 = coef_2.loc[abs_coef.index, :]

# first plot
# create dataset
y_pos = np.arange(nb_coef)[::-1]
# Create horizontal bars
plt.figure()
plt.barh(y_pos, coef_1['coef'].values)
# Create names on the x-axis
plt.yticks(y_pos, coef_1.index)
plt.vlines(x=0, ymin=0, ymax=nb_coef-1, colors='gray', linestyles='dashed')

plt.xlabel('coefficient')
plt.ylabel('cell type')
plt.title('Epidish BPRNACanProMet')
plt.savefig("../../data/processed/Deconvolution_paper_revision/all_signatures/no_bladder/l1_ratios-default/LogisticRegressionCV_coefficients_signature-EpiDISH_BPRNACanProMet_split-CV-5folds.png", bbox_inches='tight', facecolor='white')

maxi = coef_1['abs coef'].max() * 1.05
plt.xlim([-maxi, maxi])
plt.savefig("../../data/processed/Deconvolution_paper_revision/all_signatures/no_bladder/l1_ratios-default/LogisticRegressionCV_coefficients_signature-EpiDISH_BPRNACanProMet_split-CV-5folds_centered.png", bbox_inches='tight', facecolor='white')

# second plot
# create dataset
y_pos = np.arange(nb_coef)[::-1]
# Create horizontal bars
plt.figure()
plt.barh(y_pos, coef_2['coef'].values)
# Create names on the x-axis
plt.yticks(y_pos, coef_2.index)
plt.vlines(x=0, ymin=0, ymax=nb_coef-1, colors='gray', linestyles='dashed')

plt.xlabel('coefficient')
plt.ylabel('cell type')
plt.title('deconRNAseq BPRNACan')
plt.savefig("../../data/processed/Deconvolution_paper_revision/all_signatures/no_bladder/l1_ratios-default/LogisticRegressionCV_coefficients_signature-deconRNAseq_BPRNACan_split-CV-5folds.png", bbox_inches='tight', facecolor='white')

maxi = coef_2['abs coef'].max() * 1.05
plt.xlim([-maxi, maxi])
plt.savefig("../../data/processed/Deconvolution_paper_revision/all_signatures/no_bladder/l1_ratios-default/LogisticRegressionCV_coefficients_signature-deconRNAseq_BPRNACan_split-CV-5folds_centered.png", bbox_inches='tight', facecolor='white')
plt.show()

# %%
