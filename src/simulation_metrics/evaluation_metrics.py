import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score,precision_recall_curve,auc,matthews_corrcoef,f1_score

###建立灵敏度评估函数
def Sensitivity_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum):
    #result = result[result[sig_col_name]!=1]
    all_feature = result['Feature']
    result = result.dropna(subset=['p.adj.val'], axis=0,how='any')
    if result.empty:
        return 0
    if 'metadata' in result.columns:
        result = result[result['metadata']=='Group']
    predicted_feature = result['Feature'][result[sig_col_name]<pvalue]
    #predicted_feature = set(predicted_feature).intersection(set(profile.index[np.array(p_list)<0.05]))
    predicted_non_feature = set(all_feature).difference(set(predicted_feature))
    ground_truth_feature = ground_truth[(ground_truth['metadata_datum'] == metadata_datum) & (ground_truth['effect_size'].abs() >= 1)]['feature_spiked']
    true_positives = len(list(set(ground_truth_feature).intersection(set(predicted_feature))))
    false_negatives = len(set(ground_truth_feature).intersection(set(predicted_non_feature)))
    if true_positives + false_negatives == 0:
        return 0
    else:
        sensitivity = true_positives / (true_positives + false_negatives)
        return sensitivity

###建立特异性评估函数
def Specificity_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum):
    #result = result[result[sig_col_name]!=1]
    all_feature = result['Feature']
    result = result.dropna(subset=['p.adj.val'], axis=0,how='any')
    if result.empty:
        return 0
    if 'metadata' in result.columns:
        result = result[result['metadata']=='Group']
    #predicted_non_feature = result['Feature'][result['Pvalue']>=pvalue]
    predicted_feature = result['Feature'][result[sig_col_name]<pvalue]
    #predicted_feature = set(predicted_feature).intersection(set(profile.index[np.array(p_list)<0.05]))
    predicted_non_feature = set(all_feature).difference(set(predicted_feature))
    ground_truth_feature = ground_truth[(ground_truth['metadata_datum'] == metadata_datum) & (ground_truth['effect_size'].abs() >= 1)]['feature_spiked']
    ground_truth_non_feature = set(all_feature).difference(set(ground_truth_feature))
    true_negatives = len(list(set(ground_truth_non_feature).intersection(set(predicted_non_feature))))
    false_positives = len(list(set(ground_truth_non_feature).intersection(set(predicted_feature))))
    if true_negatives + false_positives == 0:
        return 0
    else:
        specificity = true_negatives / (true_negatives + false_positives)
        return specificity
###建立假阳性率评估函数
def FPR_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum):
    #result = result[result[sig_col_name]!=1]
    all_feature = result['Feature']
    result = result.dropna(subset=['p.adj.val'], axis=0,how='any')
    if result.empty:
        return 1
    if 'metadata' in result.columns:
        result = result[result['metadata']=='Group']
    #predicted_non_feature = result['Feature'][result['Pvalue']>=pvalue]
    predicted_feature = result['Feature'][result[sig_col_name]<pvalue]
    #predicted_feature = set(predicted_feature).intersection(set(profile.index[np.array(p_list)<0.05]))
    predicted_non_feature = set(all_feature).difference(set(predicted_feature))
    ground_truth_feature = ground_truth[(ground_truth['metadata_datum'] == metadata_datum) & (ground_truth['effect_size'].abs() >= 1)]['feature_spiked']
    ground_truth_non_feature = set(all_feature).difference(set(ground_truth_feature))
    false_positives = len(list(set(ground_truth_non_feature).intersection(set(predicted_feature))))
    
    FPR = false_positives / len(ground_truth_non_feature)
    return FPR

###建立准确性评估函数
def Accuracy_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum):
    #result = result[result[sig_col_name]!=1]
    all_feature = result['Feature']
    result = result.dropna(subset=['p.adj.val'], axis=0,how='any')
    if result.empty:
        return 0
    if 'metadata' in result.columns:
        result = result[result['metadata']=='Group']
    #predicted_non_feature = result['Feature'][result['Pvalue']>=pvalue]
    predicted_feature = result['Feature'][result[sig_col_name]<pvalue]
    #predicted_feature = set(predicted_feature).intersection(set(profile.index[np.array(p_list)<0.05]))
    predicted_non_feature = set(all_feature).difference(set(predicted_feature))
    ground_truth_feature = ground_truth[(ground_truth['metadata_datum'] == metadata_datum) & (ground_truth['effect_size'].abs() >= 1)]['feature_spiked']
    ground_truth_non_feature = set(all_feature).difference(set(ground_truth_feature))
    true_negatives = len(list(set(ground_truth_non_feature).intersection(set(predicted_non_feature))))
    true_positives = len(list(set(ground_truth_feature).intersection(set(predicted_feature))))
    false_positives = len(list(set(ground_truth_non_feature).intersection(set(predicted_feature))))
    false_negatives = len(set(ground_truth_feature).intersection(set(predicted_non_feature)))
    if true_negatives+true_positives+false_positives+false_negatives == 0:
        return 0
    else:
        accuracy = (true_negatives+true_positives) / (true_negatives+true_positives+false_positives+false_negatives)
        #print(true_negatives+true_positives+false_positives+false_negatives)
        return accuracy

###建立AUROC指标评估函数
def AUROC_evaluation(result,ground_truth,sig_col_name,metadata_datum):
    if 'metadata' in result.columns:
        result = result[result['metadata']=='Group']
    result.loc[:, 'p.adj.val'] = result['p.adj.val'].fillna(1)
    if result.empty:
        auc = 0
    else:
        P_value = result[sig_col_name].values
        pred_probability = 1-P_value
        #pred_probability = [pred_probability[i] if result.iloc[i,1]!=0 else 0 for i in range(len(pred_probability))]
        ground_truth_feature = ground_truth[(ground_truth['metadata_datum'] == metadata_datum) & (ground_truth['effect_size'].abs() >= 1)]['feature_spiked']
        label = [1 if i in ground_truth_feature.values.tolist() else 0 for i in result['Feature']]
        if len(np.unique(label)) < 2:
            auc = 0
        #print(P_value)
        else:
            auc = roc_auc_score(label,pred_probability)
    return auc


###建立AUPRC指标评估函数
def AUPRC_evaluation(result,ground_truth,sig_col_name,metadata_datum):
    if 'metadata' in result.columns:
        result = result[result['metadata']=='Group']
    result.loc[:, 'p.adj.val'] = result['p.adj.val'].fillna(1)
    if result.empty:
        aupr = 0
    else:
        P_value = result[sig_col_name].values
        pred_probability = 1-P_value
        #pred_probability = [pred_probability[i] if result.iloc[i,1]!=0 else 0 for i in range(len(pred_probability))]
        ground_truth_feature = ground_truth[(ground_truth['metadata_datum'] == metadata_datum) & (ground_truth['effect_size'].abs() >= 1)]['feature_spiked']
        label = [1 if i in ground_truth_feature.values.tolist() else 0 for i in result['Feature']]
        precision, recall, thresholds = precision_recall_curve(label,pred_probability)
        aupr = auc(recall,precision)
    return aupr
  



###建立precision指标函数###
def Precision_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum):
    #result = result[result[sig_col_name]!=1]
    all_feature = result['Feature']
    result = result.dropna(subset=['p.adj.val'], axis=0,how='any')
    if result.empty:
        return 0
    if 'metadata' in result.columns:
        result = result[result['metadata']=='Group']
    #predicted_non_feature = result['Feature'][result['Pvalue']>=pvalue]
    predicted_feature = result['Feature'][result[sig_col_name]<pvalue]
    #predicted_feature = set(predicted_feature).intersection(set(profile.index[np.array(p_list)<0.05]))
    predicted_non_feature = set(all_feature).difference(set(predicted_feature))
    ground_truth_feature = ground_truth[(ground_truth['metadata_datum'] == metadata_datum) & (ground_truth['effect_size'].abs() >= 1)]['feature_spiked']
    ground_truth_non_feature = set(all_feature).difference(set(ground_truth_feature))
    true_positives = len(list(set(ground_truth_feature).intersection(set(predicted_feature))))
    false_positives = len(list(set(ground_truth_non_feature).intersection(set(predicted_feature))))
    if true_positives + false_positives == 0:
        return 0
    else:
        precision = true_positives / (true_positives + false_positives)
        return precision



####建立FDR评估指标函数#####
def FDR_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum):
    #result = result[result[sig_col_name]!=1]
    all_feature = result['Feature']
    result = result.dropna(subset=['p.adj.val'], axis=0,how='any')
    if result.empty:
        return 1
    if 'metadata' in result.columns:
        result = result[result['metadata']=='Group']
    #predicted_non_feature = result['Feature'][result['Pvalue']>=pvalue]
    predicted_feature = result['Feature'][result[sig_col_name]<pvalue]
    #predicted_feature = set(predicted_feature).intersection(set(profile.index[np.array(p_list)<0.05]))
    predicted_non_feature = set(all_feature).difference(set(predicted_feature))
    ground_truth_feature = ground_truth[(ground_truth['metadata_datum'] == metadata_datum) & (ground_truth['effect_size'].abs() >= 1)]['feature_spiked']
    ground_truth_non_feature = set(all_feature).difference(set(ground_truth_feature))
    true_positives = len(list(set(ground_truth_feature).intersection(set(predicted_feature))))
    false_positives = len(list(set(ground_truth_non_feature).intersection(set(predicted_feature))))
    if true_positives + false_positives == 0:
        return 1
    else:
        FDR = false_positives / (true_positives + false_positives)
        return FDR

def F1_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum):
    #result = result[result[sig_col_name]!=1]
    all_feature = result['Feature']
    result = result.dropna(subset=['p.adj.val'], axis=0,how='any')
    if result.empty:
        return 0
    if 'metadata' in result.columns:
        result = result[result['metadata']=='Group']
    #predicted_non_feature = result['Feature'][result['Pvalue']>=pvalue]
    predicted_feature = result['Feature'][result[sig_col_name]<pvalue]
    #predicted_feature = set(predicted_feature).intersection(set(profile.index[np.array(p_list)<0.05]))
    predicted_non_feature = set(all_feature).difference(set(predicted_feature))
    ground_truth_feature = ground_truth[(ground_truth['metadata_datum'] == metadata_datum) & (ground_truth['effect_size'].abs() >= 1)]['feature_spiked']
    ground_truth_non_feature = set(all_feature).difference(set(ground_truth_feature))
    true_positives = len(list(set(ground_truth_feature).intersection(set(predicted_feature))))
    false_positives = len(list(set(ground_truth_non_feature).intersection(set(predicted_feature))))
    false_negatives = len(set(ground_truth_feature).intersection(set(predicted_non_feature)))
    if true_positives + false_negatives == 0:
        return 0
    else:
        sensitivity = true_positives / (true_positives + false_negatives)
    if true_positives + false_positives == 0:
        return 0
    else:
        precision = true_positives / (true_positives + false_positives)
    if sensitivity + precision == 0:
        return 0
    else:
        F1 = (2*sensitivity*precision)/(sensitivity+precision)
        return F1

def MCC_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum):
    if 'metadata' in result.columns:
        result = result[result['metadata']=='Group']
    result.loc[:, sig_col_name] = result[sig_col_name].fillna(1)
    if result.empty:
        return 0.0
    else:
        pred_label = (result[sig_col_name] < pvalue).astype(int).tolist()
        #pred_probability = [pred_probability[i] if result.iloc[i,1]!=0 else 0 for i in range(len(pred_probability))]
        ground_truth_feature = ground_truth[(ground_truth['metadata_datum'] == metadata_datum) & (ground_truth['effect_size'].abs() >= 1)]['feature_spiked']
        label = [1 if i in ground_truth_feature.values.tolist() else 0 for i in result['Feature']]
        mcc = matthews_corrcoef(label,pred_label)
        return mcc

def macroF1_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum):
    if 'metadata' in result.columns:
        result = result[result['metadata']=='Group']
    result.loc[:, sig_col_name] = result[sig_col_name].fillna(1)
    if result.empty:
        return 0.0
    else:
        pred_label = (result[sig_col_name] < pvalue).astype(int).tolist()
        #pred_probability = [pred_probability[i] if result.iloc[i,1]!=0 else 0 for i in range(len(pred_probability))]
        ground_truth_feature = ground_truth[(ground_truth['metadata_datum'] == metadata_datum) & (ground_truth['effect_size'].abs() >= 1)]['feature_spiked']
        label = [1 if i in ground_truth_feature.values.tolist() else 0 for i in result['Feature']]
        macrof1 = f1_score(label,pred_label,average='macro')
        return macrof1


def evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum):
    AUC = AUROC_evaluation(result,ground_truth,sig_col_name,metadata_datum)
    FPR = FPR_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum)
    Specificity = Specificity_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum)
    Sensitivity = Sensitivity_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum)
    Accuracy = Accuracy_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum)
    FDR = FDR_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum)
    Precision = Precision_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum)
    AUPR = AUPRC_evaluation(result,ground_truth,sig_col_name,metadata_datum)
    F1 = F1_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum)
    
    return FPR,Specificity,Sensitivity,Accuracy,AUC,Precision,FDR,AUPR,F1

def evaluation2(result,ground_truth,sig_col_name,pvalue,metadata_datum):
    MCC = MCC_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum)
    macroF1 = macroF1_evaluation(result,ground_truth,sig_col_name,pvalue,metadata_datum)
    return MCC,macroF1



