"""
Train ElasticNet penalized Logistic regression on deconvolution data
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
DROP_BLADDER = False
MERGE_OLD_SIG = False
ALL_SIG = True

path_interm = Path("../../data/intermediate/Deconvolution_paper")
if MERGE_OLD_SIG:
    dir_save = Path("../../data/processed/Deconvolution_paper/new_and_old_signatures")
elif ALL_SIG:
    dir_save = Path("../../data/processed/Deconvolution_paper/all_signatures")
if DROP_BLADDER:
    dir_save = dir_save / "no_bladder"
else:
    dir_save = dir_save / "all_public_datasets"
dir_save.mkdir(parents=True, exist_ok=True)

response = pd.read_csv(path_interm / "response_train.csv", index_col=0)
deconv = pd.read_csv(path_interm / "combat_all5TPM_6sign_deconvolution-cleaned.csv", index_col=0)
if MERGE_OLD_SIG:
    deconv_old = pd.read_csv(path_interm / "Gide_Maha_Riaz_Hugo_batch_corrected_deconvolution_2021-03-23-cleaned.csv", index_col=0)
    drop_col = [x for x in deconv_old.columns if (x.startswith('Quantiseq') or x.startswith('MCP'))]
    deconv_old.drop(labels=drop_col, axis=1, inplace=True)
    new_col = ['old ' + x for x in deconv_old.columns]
    deconv_old.columns = new_col
    deconv = deconv.join(deconv_old, how='inner')


common_ids = set(response.index).intersection(deconv.index)
df_all = deconv.loc[common_ids, :]
response = response.loc[df_all.index, :]

# drop variables with unique values
drop_col = [x for x in df_all.columns if df_all[x].unique().size == 1]
print("dropping these columns because they have a unique value:")
for i in drop_col:
    print("    ", i)
matches = ['XCELL', 'TAMs', 'CD226', 'Maha']
drop_col = drop_col + [x for x in df_all.columns if any(y in x for y in matches)]
if len(drop_col) > 0:
    df_all.drop(columns=drop_col, inplace=True)
# drop Maha samples
drop_row = [x for x in df_all.index if 'Maha' in x]
if DROP_BLADDER:
    drop_row = drop_row + [x for x in df_all.index if 'Snyder' in x]
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
# %% Rename columns for article

# from Ting:
# QuanTIseq and DeconRNAseq as quanTIseq and deconRNAseq 
# BPprom_314: BPRNACan3Dprom
# BPRNACan_32: BPRNAanCan3D_illumina
# BPRNACan_52: BPRNACan3Dprom_part
# BPRNACanProMet3D: BPRNACan3Dprom_enhan
# BPRNACan3DMet: BPRNACan3DMet

# we use '-' before 'illumina', 'part' and 'enhan' to latter
# avoid unintended inclusion of these variables

changes = [
    ['Quantiseq', 'quanTIseq'],   # we have Quantiseq, not QuanTIseq
    ['DeconRNA', 'deconRNAseq'],  # we have DeconRNA, not DeconRNAseq
    ['BPprom_314', 'BPRNACan3Dprom'],
    ['BPRNACan_32', 'BPRNAanCan3D-illumina'],
    ['BPRNACan_52', 'BPRNACan3Dprom-part'],
    ['BPRNACanProMet3D', 'BPRNACan3Dprom-enhan'],
    ['BPRNACan3DMet', 'BPRNACan3DMet'],
]

for c in changes:
    new_col = [x.replace(c[0], c[1]) for x in df_all.columns]
    df_all.columns = new_col

# %% Compare models trained on single signatures

if MERGE_OLD_SIG:
    conditions = [
        'Quantiseq',
        'MCP',
        'Epidish_BPprom_314',
        'DeconRNA_BPprom_314',
        'Epidish_BPRNACan',
        'DeconRNA_BPRNACan',
        'Epidish_BPRNACan3DMet',
        'DeconRNA_BPRNACan3DMet',
        'Epidish_BPRNACanProMet3D',
        'DeconRNA_BPRNACanProMet3D',
        'old Epidish_BPprom_314',
        'old DeconRNA_BPprom_314',
        'old Epidish_BPRNACan',
        'old DeconRNA_BPRNACan',
        'old Epidish_BPRNACan_32',
        'old DeconRNA_BPRNACan_32',
        'old Epidish_BPRNACan_52',
        'old DeconRNA_BPRNACan_52',
        'old Epidish_BPRNACan3DMet',
        'old DeconRNA_BPRNACan3DMet',
        'old Epidish_BPRNACanProMet3D',
        'old DeconRNA_BPRNACanProMet3D',
    ]
elif ALL_SIG:
    # we use trailing underscores to avoid including derived signatures
    # like Epidish_BPRNACan3Dprom --> Epidish_BPRNACan3Dprom-enhan
    conditions = [
        'quanTIseq',
        'MCP',
        'Epidish_BPRNACan3Dprom_',
        'deconRNAseq_BPRNACan3Dprom_',
        'Epidish_BPRNACan_',
        'deconRNAseq_BPRNACan_',
        'Epidish_BPRNACan3DMet',
        'deconRNAseq_BPRNACan3DMet',
        'Epidish_BPRNACan3Dprom-enhan',
        'deconRNAseq_BPRNACan3Dprom-enhan',
        'Epidish_BPRNAanCan3D-illumina',
        'deconRNAseq_BPRNAanCan3D-illumina',
        'Epidish_BPRNACan3Dprom-part',
        'deconRNAseq_BPRNACan3Dprom-part',
    ]
else:
    conditions = [
        'Quantiseq',
        'MCP',
        'Epidish_BPprom_314',
        'DeconRNA_BPprom_314',
        'Epidish_BPRNACan',
        'DeconRNA_BPRNACan',
        'Epidish_BPRNACan3DMet',
        'DeconRNA_BPRNACan3DMet',
        'Epidish_BPRNACanProMet3D',
        'DeconRNA_BPRNACanProMet3D',
    ]

if DROP_BLADDER or MERGE_OLD_SIG:
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
        var_idx = [x for x in df_all.columns if condi in x]
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
# %% Change QuanTIseq and DeconRNAseq to quanTIseq and deconRNAseq

coef = pd.read_csv("../../data/processed/Deconvolution_paper/all_signatures/all_public_datasets/l1_ratios-default/LogisticRegressionCV_coefficients_signature-DeconRNA_BPRNACanProMet3D_split-CV-5folds.csv", index_col=0)

new_index = [x.replace('DeconRNA', 'deconRNAseq') for x in coef.index]
coef.index = new_index

nb_coef = coef.shape[0]
nb_coef_plot = min(20, nb_coef)
labels = coef.index[:nb_coef_plot]
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
plt.savefig("../../data/processed/Deconvolution_paper/all_signatures/all_public_datasets/l1_ratios-naive/LogisticRegressionCV_coefficients_signature-deconRNAseq_BPRNACanProMet3D_split-CV-5folds.png", bbox_inches='tight', facecolor='white')
plt.show()
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

#%% Targeted comparison
condi1 = 'Epidish_BPRNACan3Dprom'
condi2 = 'deconRNAseq_BPRNACan'
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
# %% Bar plots

coef_1 = pd.read_csv("../../data/processed/Deconvolution_paper/all_signatures/all_public_datasets/l1_ratios-default/LogisticRegressionCV_coefficients_signature-Epidish_BPRNACan3Dprom_split-CV-5folds.csv", index_col=0)
new_index = [x.replace('Epidish_BPRNACan3Dprom_', '') for x in coef_1.index]
# new_index = [x.replace('BPRNACan3Dprom', 'BPRNACanProMet') for x in coef.index]
# new_index = [x.replace('DeconRNA', 'deconRNAseq') for x in coef.index]
coef_1.index = new_index
nb_coef = coef_1.shape[0]
coef_1.index.name = 'cell type'

coef_2 = pd.read_csv("../../data/processed/Deconvolution_paper/all_signatures/all_public_datasets/l1_ratios-default/LogisticRegressionCV_coefficients_signature-deconRNAseq_BPRNACan_split-CV-5folds.csv", index_col=0)
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
plt.savefig("../../data/processed/Deconvolution_paper/all_signatures/all_public_datasets/l1_ratios-default/LogisticRegressionCV_coefficients_signature-Epidish_BPRNACanProMet_split-CV-5folds.png", bbox_inches='tight', facecolor='white')

maxi = coef_1['abs coef'].max() * 1.05
plt.xlim([-maxi, maxi])
plt.savefig("../../data/processed/Deconvolution_paper/all_signatures/all_public_datasets/l1_ratios-default/LogisticRegressionCV_coefficients_signature-Epidish_BPRNACanProMet_split-CV-5folds_centered.png", bbox_inches='tight', facecolor='white')

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
plt.savefig("../../data/processed/Deconvolution_paper/all_signatures/all_public_datasets/l1_ratios-default/LogisticRegressionCV_coefficients_signature-deconRNAseq_BPRNACan_split-CV-5folds.png", bbox_inches='tight', facecolor='white')

maxi = coef_2['abs coef'].max() * 1.05
plt.xlim([-maxi, maxi])
plt.savefig("../../data/processed/Deconvolution_paper/all_signatures/all_public_datasets/l1_ratios-default/LogisticRegressionCV_coefficients_signature-deconRNAseq_BPRNACan_split-CV-5folds_centered.png", bbox_inches='tight', facecolor='white')
plt.show()
# %%
