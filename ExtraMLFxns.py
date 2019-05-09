#general requirements
import pandas as pd
import numpy as np
import sys
import itertools
import matplotlib.pyplot as plt
#process results
from sklearn.model_selection import learning_curve
from sklearn import preprocessing
from sklearn.metrics import roc_curve, auc, recall_score, f1_score, precision_score, confusion_matrix, matthews_corrcoef, hamming_loss

def subset_data(df, subset, group_split_name = None):
    X = df.drop(["Catalytic", "SITE_ID", 'cath_class', 'cath_arch', 'scop_class', 'scop_fold', 'ECOD_arch', 'ECOD_x_poshom', 'ECOD_hom'], axis = 1)
    bad_terms = ("hbond_lr_", 'dslf_fa13', 'pro_close')
    X = X.drop(columns = [term for term in X if term.startswith(bad_terms)])
    #print(X.shape, list(X))

    #general terms
    gen_set = ['MetalCodes', 'MetalAtoms', 'Depth', 'Vol']
    gen_terms = ("BSA", 'expHP', 'LoopDSSP', 'HelixDSSP', 'SheetDSSP')
    all_gen_set = [ term for term in X if term.startswith(gen_terms) ]
    gen_shell = [name for name in all_gen_set if "_S" in name]
    gen_sph = list(set(all_gen_set).difference(gen_shell))
    gen_shell += gen_set
    gen_sph += gen_set
    all_gen_set += gen_set
    #Rosetta terms only
    ros_sph1 = list(set([name for name in X if name.endswith("_3")]).difference(all_gen_set))
    ros_sph2 = list(set([ name for name in X if name.endswith("_5") ]).difference(all_gen_set))
    ros_sph3 = list(set([ name for name in X if name.endswith("_7.5") ]).difference(all_gen_set))
    ros_shell1 = list(set([ name for name in X if name.endswith("_S5") ]).difference(all_gen_set))
    ros_shell2 = list(set([ name for name in X if name.endswith("_S7.5") ]).difference(all_gen_set))
    ros_shell3 = list(set([ name for name in X if name.endswith("_S10") ]).difference(all_gen_set))
    electro = [name for name in X if name.startswith("Elec")]
    geom = [name for name in X if name.startswith("geom")]
    #pocket features only
    pocket_set = ['MetalCodes', 'MetalAtoms', 'SEPocket', 'Depth', 'Vol', 'LongPath', 'farPtLow', 'PocketAreaLow', 'OffsetLow', 'LongAxLow', 'ShortAxLow', 'farPtMid', 'PocketAreaMid', 'OffsetMid', 'LongAxMid', 'ShortAxMid', 'farPtHigh', 'PocketAreaHigh', 'OffsetHigh', 'LongAxHigh', 'ShortAxHigh']
    #pocket lining only
    lining_set = ['num_pocket_bb', 'num_pocket_sc', 'avg_eisen_hp', 'min_eisen', 'max_eisen', 'skew_eisen', 'std_dev_eisen', 'avg_kyte_hp', 'min_kyte', 'max_kyte', 'skew_kyte', 'std_dev_kyte', 'occ_vol', 'NoSC_vol', 'SC_vol_perc']

    subset_list = ["AllSph", "AllShell", "GenSph", "GenShell", "Pocket", "Lining", 
                    'RosSph', 'RosSph1','RosSph2','RosSph3', 
                    "RosShell", "RosShell1", 'RosShell2', 'RosShell3', 
                    "LinPocket", "LinRosSph", "LinRosShell", 
                    "PocketRosSph", "PocketRosShell",#new terms
                    "Geom", "LinPocketGeom", "GeomElectro", "GeomRosSph", "GeomRosShell",
                    "Electro", "ElectroRosSph", "ElectroRosShell", "LinPocketElectro", "LinPocketElectroGeom", 
                    "AllSphMinusGen", "AllSphMinusLin", "AllSphMinusPocket", 
                    "AllSphMinusGeom", "AllSphMinusElectro", "AllMinusSph", 
                    "AllShellMinusGen", "AllShellMinusLin", "AllShellMinusPocket", 
                    "AllShellMinusGeom", "AllShellMinusElectro", "AllMinusShell"
                    ]
    column_subsets = [ sorted(set(gen_sph+ros_sph1+ros_sph2+ros_sph3+pocket_set+lining_set+electro+geom)),sorted(set(gen_shell+ros_shell1+ros_shell2+ros_shell3+pocket_set+lining_set+electro+geom)), gen_sph, gen_shell, pocket_set, lining_set,
                        ros_sph1+ros_sph2+ros_sph3, ros_sph1, ros_sph2, ros_sph3,
                        ros_shell1+ros_shell2+ros_shell3, ros_shell1, ros_shell2, ros_shell3,
                        lining_set+pocket_set, lining_set+ros_sph1+ros_sph2+ros_sph3, lining_set+ros_shell1+ros_shell2+ros_shell3, 
                        pocket_set+ros_sph1+ros_sph2+ros_sph3, pocket_set+ros_shell1+ros_shell2+ros_shell3,
                        geom, lining_set+pocket_set+geom, geom+electro, geom+ros_sph1+ros_sph2+ros_sph3, geom+ros_shell1+ros_shell2+ros_shell3,
                        electro, electro+ros_sph1+ros_sph2+ros_sph3, electro+ros_shell1+ros_shell2+ros_shell3, lining_set+pocket_set+electro, lining_set+pocket_set+electro+geom, 
                        sorted(set(ros_sph1+ros_sph2+ros_sph3+pocket_set+lining_set+electro+geom)), sorted(set(gen_sph+ros_sph1+ros_sph2+ros_sph3+pocket_set+electro+geom)), sorted(set(gen_sph+ros_sph1+ros_sph2+ros_sph3+lining_set+electro+geom)),
                        sorted(set(gen_sph+ros_sph1+ros_sph2+ros_sph3+pocket_set+lining_set+electro)), sorted(set(gen_sph+ros_sph1+ros_sph2+ros_sph3+pocket_set+lining_set+geom)), sorted(set(gen_sph+pocket_set+lining_set+electro+geom)),
                        sorted(set(ros_shell1+ros_shell2+ros_shell3+pocket_set+lining_set+electro+geom)), sorted(set(gen_shell+ros_shell1+ros_shell2+ros_shell3+pocket_set+electro+geom)), sorted(set(gen_shell+ros_shell1+ros_shell2+ros_shell3+lining_set+electro+geom)),
                        sorted(set(gen_shell+ros_shell1+ros_shell2+ros_shell3+pocket_set+lining_set+electro)), sorted(set(gen_shell+ros_shell1+ros_shell2+ros_shell3+pocket_set+lining_set+geom)), sorted(set(gen_shell+pocket_set+lining_set+electro+geom))
                        ]

    #print(column_subsets[subset_list.index(data_subset)] )
    if subset in subset_list:
        X = X[ column_subsets[subset_list.index(subset)] ] 
    else:
        print("Not a subset in list; defaulting to AllSph")
        X = X[ column_subsets[0] ] #this is all for usage with PCA/UMAP; it uses the rosetta sphere terms plus all the non-rosetta terms
    if group_split_name != None:
        X["groupID"] = preprocessing.LabelEncoder().fit_transform( df[[group_split_name]].astype(str) ) #add fold identifiers converted to number
    y = df[["Catalytic", "SITE_ID"]] #keep SITE_ID in the answers for merging back later
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
