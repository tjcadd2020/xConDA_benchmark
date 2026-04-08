import pandas as pd
from evaluation_metrics import evaluation,evaluation2
import os
import numpy as np
from sklearn.metrics import roc_auc_score,precision_recall_curve,auc
import re

#base_directory = '/home/data/ZXY/meta_causal/benchmark_simulated_data/'

#ground_truth = pd.read_csv('~/meta_causal/simulated_data/1/ground_truth1.tsv',sep='\t')

def run_evaluation(significance_threshold,metadata_datum,base_directory,num_simulated_data_start,num_simulated_data_end):
    for i in range(num_simulated_data_start,num_simulated_data_end):
        result_folder_name = str(i) + '/result/'
        folder_name = str(i)
        result_folder_path = os.path.join(base_directory,result_folder_name)
        folder_path = os.path.join(base_directory,folder_name)
        #获取相应的文件
        files = os.listdir(result_folder_path)
        folder_results = []
        
        for filename in files:
            method_result_file_path = os.path.join(result_folder_path,filename)
            method_result = pd.read_csv(method_result_file_path,sep='\t')
            match = re.search(r'rep(\d+)', filename)
            j = int(match.group(1))
            ground_truth_file_name = 'ground_truth' + str(i) + '_rep'+ str(j) + '.tsv'
            ground_truth_path = os.path.join(folder_path,ground_truth_file_name)
            ground_truth = pd.read_csv(ground_truth_path,sep='\t')
            print(filename)
            method_performance = evaluation(method_result,ground_truth,'p.adj.val',significance_threshold,metadata_datum)
            folder_results.append((filename,method_performance))
        evaluation_result_filename = 'metrics_evaluation_result_' + str(significance_threshold) + '.tsv'
        evaluation_file_path = os.path.join(folder_path,evaluation_result_filename)
        with open(evaluation_file_path,'w',encoding = 'utf-8') as eva_file:
            eva_file.write('\t'.join(['Method_rep','FPR','Specificity','Sensitivity','Accuracy','AUC','Precision','FDR','AUPR','F1']) + '\n') 
            for filename, method_performance in folder_results:
                eva_file.write(filename + '\t' +'\t'.join(map(str,method_performance)) + '\n')


####显著性P值的列名被更改为p.adj.val

def run_evaluation_ncc(significance_threshold,metadata_datum,base_directory,num_simulated_data_start,num_simulated_data_end):
    for i in range(num_simulated_data_start,num_simulated_data_end):
        result_folder_name = str(i) + '/ncc_result/'
        folder_name = str(i)
        result_folder_path = os.path.join(base_directory,result_folder_name)
        folder_path = os.path.join(base_directory,folder_name)
        #获取相应的文件
        files = os.listdir(result_folder_path)
        folder_results = []

        for filename in files:
            method_result_file_path = os.path.join(result_folder_path,filename)
            method_result = pd.read_csv(method_result_file_path,sep='\t')
            match = re.search(r'rep(\d+)', filename)
            j = int(match.group(1))
            ground_truth_file_name = 'ground_truth' + str(i) + '_rep'+ str(j) + '.tsv'
            ground_truth_path = os.path.join(folder_path,ground_truth_file_name)
            ground_truth = pd.read_csv(ground_truth_path,sep='\t')
            print(filename)
            method_performance = evaluation(method_result,ground_truth,'p.adj.val',significance_threshold,metadata_datum)
            folder_results.append((filename,method_performance))
        evaluation_result_filename = 'metrics_evaluation_result_ncc_' + str(significance_threshold) + '.tsv'
        evaluation_file_path = os.path.join(folder_path,evaluation_result_filename)
        with open(evaluation_file_path,'w',encoding = 'utf-8') as eva_file:
            eva_file.write('\t'.join(['Method_rep','FPR','Specificity','Sensitivity','Accuracy','AUC','Precision','FDR','AUPR','F1']) + '\n') 
            for filename, method_performance in folder_results:
                eva_file.write(filename + '\t' +'\t'.join(map(str,method_performance)) + '\n')

def run_evaluation_prevalence(significance_threshold,metadata_datum,base_directory,num_simulated_data_start,num_simulated_data_end,feature_prevalence_min,feature_prevalence_max):
    for i in range(num_simulated_data_start,num_simulated_data_end):
        result_folder_name = str(i) + '/result/'
        folder_name = str(i)
        result_folder_path = os.path.join(base_directory,result_folder_name)
        folder_path = os.path.join(base_directory,folder_name)
        #获取相应的文件
        files = os.listdir(result_folder_path)
        folder_results = []

        for filename in files:
            method_result_file_path = os.path.join(result_folder_path,filename)
            method_result = pd.read_csv(method_result_file_path,sep='\t')
            match = re.search(r'rep(\d+)', filename)
            j = int(match.group(1))
            simulated_data_name = 'simulated_data_abs_' + str(i) + '_rep' + str(j) +'.tsv'
            simulated_data_path =  os.path.join(folder_path,simulated_data_name)
            simulated_data = pd.read_csv(simulated_data_path,sep='\t')
            prevalence = (simulated_data != 0).sum(axis=0) / len(simulated_data)
            prevalence_feature = prevalence[(prevalence > feature_prevalence_min) & (prevalence <= feature_prevalence_max)].index
            method_result = method_result[method_result["Feature"].isin(prevalence_feature)]

            ground_truth_file_name = 'ground_truth' + str(i) + '_rep'+ str(j) + '.tsv'
            ground_truth_path = os.path.join(folder_path,ground_truth_file_name)
            ground_truth = pd.read_csv(ground_truth_path,sep='\t')
            ground_truth = ground_truth[ground_truth["feature_spiked"].isin(prevalence_feature)]
            print(filename)
            method_performance = evaluation(method_result,ground_truth,'p.adj.val',significance_threshold,metadata_datum)
            folder_results.append((filename,method_performance))
        evaluation_result_filename = 'metrics_evaluation_result_prevalence' + str(feature_prevalence_min) + '-' + str(feature_prevalence_max) + '_' + str(significance_threshold) + '.tsv'
        evaluation_file_path = os.path.join(folder_path,evaluation_result_filename)
        with open(evaluation_file_path,'w',encoding = 'utf-8') as eva_file:
            eva_file.write('\t'.join(['Method_rep','FPR','Specificity','Sensitivity','Accuracy','AUC','Precision','FDR','AUPR','F1']) + '\n') 
            for filename, method_performance in folder_results:
                eva_file.write(filename + '\t' +'\t'.join(map(str,method_performance)) + '\n')

def run_evaluation_prevalence2(significance_threshold,metadata_datum,base_directory,num_simulated_data_start,num_simulated_data_end,feature_prevalence_min,feature_prevalence_max):
    for i in range(num_simulated_data_start,num_simulated_data_end):
        result_folder_name = str(i) + '/result/'
        folder_name = str(i)
        result_folder_path = os.path.join(base_directory,result_folder_name)
        folder_path = os.path.join(base_directory,folder_name)
        #获取相应的文件
        files = os.listdir(result_folder_path)
        folder_results = []

        for filename in files:
            method_result_file_path = os.path.join(result_folder_path,filename)
            method_result = pd.read_csv(method_result_file_path,sep='\t')
            match = re.search(r'rep(\d+)', filename)
            j = int(match.group(1))
            simulated_data_name = 'simulated_data_abs_' + str(i) + '_rep' + str(j) +'.tsv'
            simulated_data_path =  os.path.join(folder_path,simulated_data_name)
            simulated_data = pd.read_csv(simulated_data_path,sep='\t')
            prevalence = (simulated_data != 0).sum(axis=0) / len(simulated_data)
            prevalence_feature = prevalence[(prevalence > feature_prevalence_min) & (prevalence <= feature_prevalence_max)].index
            method_result = method_result[method_result["Feature"].isin(prevalence_feature)]

            ground_truth_file_name = 'ground_truth' + str(i) + '_rep'+ str(j) + '.tsv'
            ground_truth_path = os.path.join(folder_path,ground_truth_file_name)
            ground_truth = pd.read_csv(ground_truth_path,sep='\t')
            ground_truth = ground_truth[ground_truth["feature_spiked"].isin(prevalence_feature)]
            print(filename)
            method_performance = evaluation2(method_result,ground_truth,'p.adj.val',significance_threshold,metadata_datum)
            folder_results.append((filename,method_performance))
        evaluation_result_filename = 'MCC_evaluation_result_prevalence' + str(feature_prevalence_min) + '-' + str(feature_prevalence_max) + '_' + str(significance_threshold) + '.tsv'
        evaluation_file_path = os.path.join(folder_path,evaluation_result_filename)
        with open(evaluation_file_path,'w',encoding = 'utf-8') as eva_file:
            eva_file.write('\t'.join(['Method_rep','MCC','Macro-F1']) + '\n') 
            for filename, method_performance in folder_results:
                eva_file.write(filename + '\t' +'\t'.join(map(str,method_performance)) + '\n')

def run_evaluation2(significance_threshold,metadata_datum,base_directory,num_simulated_data_start,num_simulated_data_end):
    for i in range(num_simulated_data_start,num_simulated_data_end):
        result_folder_name = str(i) + '/result/'
        folder_name = str(i)
        result_folder_path = os.path.join(base_directory,result_folder_name)
        folder_path = os.path.join(base_directory,folder_name)
        #获取相应的文件
        files = os.listdir(result_folder_path)
        folder_results = []
        
        for filename in files:
            method_result_file_path = os.path.join(result_folder_path,filename)
            method_result = pd.read_csv(method_result_file_path,sep='\t')
            match = re.search(r'rep(\d+)', filename)
            j = int(match.group(1))
            ground_truth_file_name = 'ground_truth' + str(i) + '_rep'+ str(j) + '.tsv'
            ground_truth_path = os.path.join(folder_path,ground_truth_file_name)
            ground_truth = pd.read_csv(ground_truth_path,sep='\t')
            print(filename)
            method_performance = evaluation2(method_result,ground_truth,'p.adj.val',significance_threshold,metadata_datum)
            folder_results.append((filename,method_performance))
        evaluation_result_filename = 'MCC_evaluation_result_' + str(significance_threshold) + '.tsv'
        evaluation_file_path = os.path.join(folder_path,evaluation_result_filename)
        with open(evaluation_file_path,'w',encoding = 'utf-8') as eva_file:
            eva_file.write('\t'.join(['Method_rep','MCC','Macro-F1']) + '\n') 
            for filename, method_performance in folder_results:
                #print(method_performance)
                eva_file.write(filename + '\t' +'\t'.join(map(str,method_performance)) + '\n')


def run_evaluation_meta(significance_threshold,metadata_datum,base_directory,num_simulated_data_start,num_simulated_data_end):
    for i in range(num_simulated_data_start,num_simulated_data_end):
        result_folder_name = str(i) + '/result/'
        folder_name = str(i)
        result_folder_path = os.path.join(base_directory,result_folder_name)
        folder_path = os.path.join(base_directory,folder_name)
        #获取相应的文件
        files = os.listdir(result_folder_path)
        folder_results = []
        
        for filename in files:
            method_result_file_path = os.path.join(result_folder_path,filename)
            method_result = pd.read_csv(method_result_file_path,sep='\t')
            match = re.search(r'rep(\d+)', filename)
            j = int(match.group(1))
            ground_truth_file_name = 'ground_truthmeta' + '_rep'+ str(j) + '.tsv'
            ground_truth_path = os.path.join(folder_path,ground_truth_file_name)
            ground_truth = pd.read_csv(ground_truth_path,sep='\t')
            print(filename)
            method_performance = evaluation(method_result,ground_truth,'p.adj.val',significance_threshold,metadata_datum)
            folder_results.append((filename,method_performance))
        evaluation_result_filename = 'metrics_evaluation_result_' + str(significance_threshold) + '.tsv'
        evaluation_file_path = os.path.join(folder_path,evaluation_result_filename)
        with open(evaluation_file_path,'w',encoding = 'utf-8') as eva_file:
            eva_file.write('\t'.join(['Method_rep','FPR','Specificity','Sensitivity','Accuracy','AUC','Precision','FDR','AUPR','F1']) + '\n') 
            for filename, method_performance in folder_results:
                eva_file.write(filename + '\t' +'\t'.join(map(str,method_performance)) + '\n')

    

def run_evaluation_meta2(significance_threshold,metadata_datum,base_directory,num_simulated_data_start,num_simulated_data_end):
    for i in range(num_simulated_data_start,num_simulated_data_end):
        result_folder_name = str(i) + '/result/'
        folder_name = str(i)
        result_folder_path = os.path.join(base_directory,result_folder_name)
        folder_path = os.path.join(base_directory,folder_name)
        #获取相应的文件
        files = os.listdir(result_folder_path)
        folder_results = []
        
        for filename in files:
            method_result_file_path = os.path.join(result_folder_path,filename)
            method_result = pd.read_csv(method_result_file_path,sep='\t')
            match = re.search(r'rep(\d+)', filename)
            j = int(match.group(1))
            ground_truth_file_name = 'ground_truthmeta' + '_rep'+ str(j) + '.tsv'
            ground_truth_path = os.path.join(folder_path,ground_truth_file_name)
            ground_truth = pd.read_csv(ground_truth_path,sep='\t')
            print(filename)
            method_performance = evaluation2(method_result,ground_truth,'p.adj.val',significance_threshold,metadata_datum)
            folder_results.append((filename,method_performance))
        evaluation_result_filename = 'MCC_evaluation_result_' + str(significance_threshold) + '.tsv'
        evaluation_file_path = os.path.join(folder_path,evaluation_result_filename)
        with open(evaluation_file_path,'w',encoding = 'utf-8') as eva_file:
            eva_file.write('\t'.join(['Method_rep','MCC','Macro-F1']) + '\n') 
            for filename, method_performance in folder_results:
                #print(method_performance)
                eva_file.write(filename + '\t' +'\t'.join(map(str,method_performance)) + '\n')


