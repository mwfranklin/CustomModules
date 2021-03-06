#general requirements
import pandas as pd
import numpy as np
import sys
import itertools
from collections import Counter
import matplotlib.pyplot as plt
#process results
from sklearn.model_selection import learning_curve
from sklearn import preprocessing
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_curve, auc, recall_score, f1_score, precision_score, confusion_matrix, matthews_corrcoef, hamming_loss, jaccard_score

def subset_data(df, subset, group_split_name = None, other_bad_terms = []):
    not_needed = ("Catalytic", "SITE_ID", "ValidSet", 'NewSet', 'cath_class', 'cath_arch', 'scop_class', 'scop_fold', 'ECOD_arch', 'ECOD_x_poshom', 'ECOD_hom')
    X = df.drop(columns = [term for term in df if term.startswith(not_needed)])
    bad_terms = ("hbond_lr_", 'dslf_fa13', 'pro_close', 'ref', 'fa_sol_', 'MetalCodes', 'MetalAtoms', 'SEPocket')
    X = X.drop(columns = [term for term in X if term.startswith(bad_terms)])
    #print(X.shape)#, list(X))
    #general terms
    gen_set = ['Depth', 'Vol', "SITEDistCenter", "SITEDistNormCenter"]
    gen_terms = ("BSA", 'SASA')
    all_gen_set = [ term for term in X if term.startswith(gen_terms) ]
    gen_shell = [name for name in all_gen_set if "_S" in name]
    gen_sph = list(set(all_gen_set).difference(gen_shell))
    gen_shell += gen_set
    gen_shell += ["BSA_3.5", "SASA_3.5"]
    gen_sph += gen_set
    all_gen_set += gen_set
    all_gen_set = sorted(set(all_gen_set))
    #Rosetta terms only
    ros_sum_sph0 = list(set([name for name in X if name.endswith("_Sum_3.5")]).difference(all_gen_set))
    ros_sum_sph1 = list(set([ name for name in X if name.endswith("_Sum_5") ]).difference(all_gen_set))
    ros_sum_sph2 = list(set([ name for name in X if name.endswith("_Sum_7.5") ]).difference(all_gen_set))
    ros_sum_sph3 = list(set([ name for name in X if name.endswith("_Sum_9") ]).difference(all_gen_set))
    ros_sum_shell1 = list(set([ name for name in X if name.endswith("_Sum_S5") ]).difference(all_gen_set))
    ros_sum_shell2 = list(set([ name for name in X if name.endswith("_Sum_S7.5") ]).difference(all_gen_set))
    ros_sum_shell3 = list(set([ name for name in X if name.endswith("_Sum_S9") ]).difference(all_gen_set))
    ros_sum_shell = ros_sum_sph0 + ros_sum_shell1 + ros_sum_shell2 + ros_sum_shell3
    ros_sum_sph = ros_sum_sph0 + ros_sum_sph1 + ros_sum_sph2 + ros_sum_sph3
    
    ros_mean_sph0 = list(set([name for name in X if name.endswith("_Mean_3.5")]).difference(all_gen_set))
    ros_mean_sph1 = list(set([ name for name in X if name.endswith("_Mean_5") ]).difference(all_gen_set))
    ros_mean_sph2 = list(set([ name for name in X if name.endswith("_Mean_7.5") ]).difference(all_gen_set))
    ros_mean_sph3 = list(set([ name for name in X if name.endswith("_Mean_9") ]).difference(all_gen_set))
    ros_mean_shell1 = list(set([ name for name in X if name.endswith("_Mean_S5") ]).difference(all_gen_set))
    ros_mean_shell2 = list(set([ name for name in X if name.endswith("_Mean_S7.5") ]).difference(all_gen_set))
    ros_mean_shell3 = list(set([ name for name in X if name.endswith("_Mean_S9") ]).difference(all_gen_set))
    ros_mean_shell = ros_mean_sph0 + ros_mean_shell1 + ros_mean_shell2 + ros_mean_shell3
    ros_mean_sph = ros_mean_sph0 + ros_mean_sph1 + ros_mean_sph2 + ros_mean_sph3
    
    electro = [name for name in X if name.startswith("Elec")]
    geom = [name for name in X if name.startswith("geom")]
    findgeo_geoms = ("lin", "trv", "tri", "tev", "spv", 
        "tet", "spl", "bva", "bvp", "pyv", 
        "spy", "tbp", "tpv", 
        "oct", "tpr", "pva", "pvp", "cof", "con", "ctf", "ctn",
        "pbp", "coc", "ctp", "hva", "hvp", "cuv", "sav",
        "hbp", "cub", "sqa", "boc", "bts", "btt", 
        "ttp", "csa")
    geom = [name for name in geom if not name.endswith(findgeo_geoms)] #remove the individual geom types
    #pocket features only
    pocket_set = ['Depth', 'Vol', "SITEDistCenter", "SITEDistNormCenter", 'LongPath', 'farPtLow', 'PocketAreaLow', 'OffsetLow', 'LongAxLow', 'ShortAxLow', 'farPtMid', 'PocketAreaMid', 'OffsetMid', 'LongAxMid', 'ShortAxMid', 'farPtHigh', 'PocketAreaHigh', 'OffsetHigh', 'LongAxHigh', 'ShortAxHigh']
    pocket_set = list(set(pocket_set).difference(other_bad_terms))
    #pocket lining only
    lining_set = ['num_pocket_bb', 'num_pocket_sc', 'avg_eisen_hp', 'min_eisen', 'max_eisen', 'skew_eisen', 'std_dev_eisen', 'avg_kyte_hp', 'min_kyte', 'max_kyte', 'skew_kyte', 'std_dev_kyte', 'occ_vol', 'NoSC_vol', 'SC_vol_perc', 'LiningArea']
    lining_set = list(set(lining_set).difference(other_bad_terms))
    
    #print( len(all_gen_set), len(sorted(set(ros_sum_sph+ros_sum_shell+ros_mean_shell+ros_mean_sph))), len(electro), len(geom), len(pocket_set), len(lining_set))
    #print(len(sorted(set(ros_sum_sph+ros_sum_shell+ros_mean_shell+ros_mean_sph+all_gen_set+electro+geom+pocket_set+lining_set))))
    
    subset_list = ["AllSumSph", "AllMeanSph", "AllSumShell", "AllMeanShell", 
                    "GenSph", "GenShell", "Pocket", "Lining", 
                    'RosSumSph', 'RosSumSph0', 'RosSumSph1', 'RosMeanSph', 'RosMeanSph0', 'RosMeanSph1', "RosSumSphInner2", "RosMeanSphInner2",
                    'RosSumShell', 'RosSumShell1', 'RosMeanShell', 'RosMeanShell1',"RosSumShellInner2", "RosMeanShellInner2",
                    "LinPocket", "LinRosSumSph", "LinRosMeanSph", "LinRosSumShell", "LinRosMeanShell",
                    "PocketRosSumSph", "PocketRosMeanSph", "PocketRosSumShell", "PocketRosMeanShell",
                    "Geom", "LinPocketGeom", "GeomElectro", "GeomRosSumSph", "GeomRosSumShell", "GeomRosMeanSph", "GeomRosMeanShell",
                    
                    "Electro", "LinPocketElectro", "LinPocketElectroGeom", "ElectroRosSumSph", "ElectroRosSumShell", "ElectroRosMeanSph", "ElectroRosMeanShell", 
                    "AllSumSphMinusGen", "AllSumSphMinusLin", "AllSumSphMinusPocket", 
                    "AllSumSphMinusGeom", "AllSumSphMinusElectro", 
                    "AllMeanSphMinusGen", "AllMeanSphMinusLin", "AllMeanSphMinusPocket", 
                    "AllMeanSphMinusGeom", "AllMeanSphMinusElectro", "AllMinusRosSph",
                    
                    "AllSumShellMinusGen", "AllSumShellMinusLin", "AllSumShellMinusPocket", 
                    "AllSumShellMinusGeom", "AllSumShellMinusElectro", 
                    "AllMeanShellMinusGen", "AllMeanShellMinusLin", "AllMeanShellMinusPocket", 
                    "AllMeanShellMinusGeom", "AllMeanShellMinusElectro", "AllMinusRosShell",
                    ]
    column_subsets = [ sorted(set(gen_sph+ros_sum_sph+pocket_set+lining_set+electro+geom)),sorted(set(gen_shell+ros_mean_sph+pocket_set+lining_set+electro+geom)), sorted(set(gen_sph+ros_sum_shell+pocket_set+lining_set+electro+geom)),sorted(set(gen_shell+ros_mean_shell+pocket_set+lining_set+electro+geom)), 
                        gen_sph, gen_shell, pocket_set, lining_set,
                        ros_sum_sph, ros_sum_sph0, ros_sum_sph1, ros_mean_sph, ros_mean_sph0, ros_mean_sph1, sorted(set(ros_sum_sph0+ros_sum_sph1)), sorted(set(ros_mean_sph0+ros_mean_sph1)),
                        ros_sum_shell, ros_sum_shell1, ros_mean_shell, ros_mean_shell1, sorted(set(ros_sum_sph0 + ros_sum_shell1)), sorted(set(ros_mean_sph0 + ros_mean_shell1)), 
                        lining_set+pocket_set, lining_set+ros_sum_sph, lining_set+ros_mean_sph, lining_set+ros_sum_shell, lining_set+ros_mean_shell, 
                        pocket_set+ros_sum_sph, pocket_set+ros_mean_sph, pocket_set+ros_sum_shell, pocket_set+ros_mean_shell, 
                        geom, lining_set+pocket_set+geom, geom+electro, geom+ros_sum_sph, geom+ros_sum_shell,geom+ros_mean_sph, geom+ros_mean_shell,
                        electro, lining_set+pocket_set+electro, lining_set+pocket_set+electro+geom, electro+ros_sum_sph, electro+ros_sum_shell, electro+ros_mean_sph, electro+ros_mean_shell, 
                        
                        sorted(set(ros_sum_sph+pocket_set+lining_set+electro+geom)), sorted(set(gen_sph+ros_sum_sph+pocket_set+electro+geom)), sorted(set(gen_sph+ros_sum_sph+lining_set+electro+geom)),
                        sorted(set(gen_sph+ros_sum_sph+pocket_set+lining_set+electro)), sorted(set(gen_sph+ros_sum_sph+pocket_set+lining_set+geom)), 
                        sorted(set(ros_mean_sph+pocket_set+lining_set+electro+geom)), sorted(set(gen_sph+ros_mean_sph+pocket_set+electro+geom)), sorted(set(gen_sph+ros_mean_sph+lining_set+electro+geom)),
                        sorted(set(gen_sph+ros_mean_sph+pocket_set+lining_set+electro)), sorted(set(gen_sph+ros_mean_sph+pocket_set+lining_set+geom)), sorted(set(gen_sph+pocket_set+lining_set+electro+geom)),
                        
                        sorted(set(ros_sum_shell+pocket_set+lining_set+electro+geom)), sorted(set(gen_shell+ros_sum_shell+pocket_set+electro+geom)), sorted(set(gen_shell+ros_sum_shell+lining_set+electro+geom)),
                        sorted(set(gen_shell+ros_sum_shell+pocket_set+lining_set+electro)), sorted(set(gen_shell+ros_sum_shell+pocket_set+lining_set+geom)), 
                        sorted(set(ros_mean_shell+pocket_set+lining_set+electro+geom)), sorted(set(gen_shell+ros_mean_shell+pocket_set+electro+geom)), sorted(set(gen_shell+ros_mean_shell+lining_set+electro+geom)),
                        sorted(set(gen_shell+ros_mean_shell+pocket_set+lining_set+electro)), sorted(set(gen_shell+ros_mean_shell+pocket_set+lining_set+geom)), sorted(set(gen_shell+pocket_set+lining_set+electro+geom)),
                        ]

    #print(column_subsets[subset_list.index(data_subset)] )
    if subset in subset_list:
        X = X[ column_subsets[subset_list.index(subset)] ] 
    else:
        print("Not a subset in list; defaulting to AllSph")
        X = X[ column_subsets[0] ] #this is all for usage with PCA/UMAP; it uses the rosetta sphere terms plus all the non-rosetta terms
    if group_split_name != None:
        X["groupID"] = preprocessing.LabelEncoder().fit_transform( df[[group_split_name]].astype(str) ) #add fold identifiers converted to number

    y = df[["SITE_ID", "Catalytic"]] #keep SITE_ID in the answers for merging back later
    return(X, y)

def plot_confusion_matrix(cm, classes, normalize=False, title='Confusion matrix', filename = "test.png", cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt),
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")

    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    plt.tight_layout()
    plt.savefig(filename)

def plot_learning_curve(estimator, title, X, y, ylim=None, cv=None, filename = "test.png", n_jobs=None, train_sizes=np.linspace(.1, 1.0, 5)):
    """
    Generate a simple plot of the test and training learning curve.

    Parameters
    ----------
    estimator : object type that implements the "fit" and "predict" methods
        An object of that type which is cloned for each validation.

    title : string
        Title for the chart.

    X : array-like, shape (n_samples, n_features)
        Training vector, where n_samples is the number of samples and
        n_features is the number of features.

    y : array-like, shape (n_samples) or (n_samples, n_features), optional
        Target relative to X for classification or regression;
        None for unsupervised learning.

    ylim : tuple, shape (ymin, ymax), optional
        Defines minimum and maximum yvalues plotted.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:
          - None, to use the default 3-fold cross-validation,
          - integer, to specify the number of folds.
          - :term:`CV splitter`,
          - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, if ``y`` is binary or multiclass,
        :class:`StratifiedKFold` used. If the estimator is not a classifier
        or if ``y`` is neither binary nor multiclass, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validators that can be used here.

    n_jobs : int or None, optional (default=None)
        Number of jobs to run in parallel.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    train_sizes : array-like, shape (n_ticks,), dtype float or int
        Relative or absolute numbers of training examples that will be used to
        generate the learning curve. If the dtype is float, it is regarded as a
        fraction of the maximum size of the training set (that is determined
        by the selected validation method), i.e. it has to be within (0, 1].
        Otherwise it is interpreted as absolute sizes of the training sets.
        Note that for classification the number of samples usually have to
        be big enough to contain at least one sample from each class.
        (default: np.linspace(0.1, 1.0, 5))
    """
    plt.figure()
    plt.title(title)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel("Training examples")
    plt.ylabel("Score")
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X, y, cv=cv, n_jobs=n_jobs, train_sizes=train_sizes)
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)
    plt.grid()

    plt.fill_between(train_sizes, train_scores_mean - train_scores_std,
                     train_scores_mean + train_scores_std, alpha=0.1,
                     color="r")
    plt.fill_between(train_sizes, test_scores_mean - test_scores_std,
                     test_scores_mean + test_scores_std, alpha=0.1, color="g")
    plt.plot(train_sizes, train_scores_mean, 'o-', color="r",
             label="Training score")
    plt.plot(train_sizes, test_scores_mean, 'o-', color="g",
             label="Cross-validation score")

    plt.legend(loc="best")
    plt.savefig(filename)
    return plt

def plot_roc_curve(names, fpr, tpr, filename):
    fig, ax1 = plt.subplots(nrows = 1, ncols= 1, figsize = (10,10))
    lw = 2
    color_list = ['darkorange','magenta','limegreen','darkorchid','cyan','deepskyblue','darkred','darkgoldenrod', 'darkgreen','grey','rosybrown', 'mediumslateblue', 'indigo', 'lightpink']
    for name in names:
        #print(name)
        ax1.plot(fpr[name], tpr[name], color=color_list[names.index(name)],lw=lw, label='%s (area = %0.2f)' %(name, roc_auc[name]))
    ax1.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    ax1.set_xlim([0.0, 1.0])
    ax1.set_ylim([0.0, 1.05])
    ax1.set_xlabel('False Positive Rate')
    ax1.set_ylabel('True Positive Rate')
    plt.title('ROC Curve')
    plt.legend(loc="lower right")
    plt.savefig(filename)

def plot_roc_curve_pickle(names, dict_values, scores, filename):
    #fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols= 2, figsize = (20,10))
    fig, ax1 = plt.subplots(nrows = 1, ncols= 1, figsize = (10,10))
    lw = 2
    color_list = ['darkorange','magenta','limegreen','darkorchid','cyan','deepskyblue','darkred','darkgoldenrod', 'darkgreen','grey','rosybrown', 'mediumslateblue', 'indigo', 'lightpink']
    for name in names:
        #print(name)
        #print(scores["AUC"][scores["Classifier"] == name] )
        ax1.plot(dict_values[name][0], dict_values[name][1], color=color_list[names.index(name)],lw=lw, label='%s (area = %0.2f)' %(name, scores["AUC"][scores["Classifier"] == name] ))
        #ax2.plot(dict_values[name][2], dict_values[name][3], color=color_list[names.index(name)],lw=lw, label='%s (area = %0.2f)' %(name, scores["AUCNeg"][scores["Classifier"] == name]))
                
    ax1.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    ax1.set_xlim([0.0, 1.0])
    ax1.set_ylim([0.0, 1.05])
    ax1.set_xlabel('False Positive Rate')
    ax1.set_ylabel('True Positive Rate')
    ax1.legend()
    """ax2.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    ax2.set_xlim([0.0, 1.0])
    ax2.set_ylim([0.0, 1.05])
    ax2.set_xlabel('False Negative Rate')
    ax2.set_ylabel('True Negative Rate')
    ax2.legend(loc="lower right")"""  
    
    plt.title('ROC Curve')
    plt.savefig(filename)

def lr_plus(y_true, y_pred):
    cnf_matrix = confusion_matrix(y_true, y_pred)
    cnf_matrix = cnf_matrix.astype('float') / cnf_matrix.sum(axis=1)[:, np.newaxis]
    #print(cnf_matrix)
    try:
        return(cnf_matrix[0][0]/ (1-cnf_matrix[1][1]) )
    except IndexError:
        return(0)

def tnr_score(y_true, y_pred):
    cnf_matrix = confusion_matrix(y_true, y_pred)
    cnf_matrix = cnf_matrix.astype('float') / cnf_matrix.sum(axis=1)[:, np.newaxis]
    #print(cnf_matrix)
    try:
        return(cnf_matrix[1][1])
    except IndexError:
        return(0)
    
def roc_area(y_true, y_pred, this_label = 0):
    fpr, tpr, _ = roc_curve(y_true, y_pred, pos_label = this_label)
    roc_auc = auc(fpr, tpr)
    return(roc_auc)

def prec_score_custom(y_true, y_pred, this_label = 0):
    return( precision_score(y_true, y_pred, pos_label= this_label) )

def mcc_score(y_true, y_pred):
    return( matthews_corrcoef(y_true, y_pred))
    
def jac_score(y_true, y_pred, this_label = 0):
    return( jaccard_score(y_true, y_pred, pos_label=this_label))    

def hamming_score(y_true, y_pred):
    return( hamming_loss(y_true, y_pred))

def dist_rand(y_true, y_pred):
    cnf_matrix = confusion_matrix(y_true, y_pred)
    totals = cnf_matrix.sum(axis=0)[:, np.newaxis]
    totals = totals/np.sum(totals)
    cnf_matrix = cnf_matrix.astype('float') / cnf_matrix.sum(axis=1)[:, np.newaxis]
    dist_rand = (cnf_matrix[0][0] + -1*(1- cnf_matrix[1][1]) )/np.sqrt(2)
    return(dist_rand)

def calc_bal_acc(y_true, y_pred):
    cnf_matrix = confusion_matrix(y_true, y_pred)
    totals = cnf_matrix.sum(axis=0)[:, np.newaxis]
    totals = totals/np.sum(totals)
    #print(totals, "\n", cnf_matrix)
    cnf_matrix = cnf_matrix.astype('float') / cnf_matrix.sum(axis=1)[:, np.newaxis]
    bal_acc = cnf_matrix[0][0] * totals[0] + cnf_matrix[1][1] * totals[1]
    #print(bal_acc)
    return(bal_acc[0])
