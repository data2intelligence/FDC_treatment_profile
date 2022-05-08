#!/usr/bin/env python

import os, sys, re, numpy, copy, pandas, requests, browsercookie
from glob import glob

# please change this to other browsers if you are not using chrome
cj = browsercookie.chrome()


def download_file_list(f, out_path):
    """ Download all sample annotations and expression data matrices (if available from FDC automatic processing from public repositories) """
    fin = open(f)
    for l in fin:
        url = l.strip()
        
        fields = url.split('/')
        target = fields[-1]
        
        # also add curator ID for meta annotations
        if target.find('.meta') > 0: target += '.' + fields[-2]
        
        print(target)

        fout = open(os.path.join(out_path, target), 'wb')
        
        read = requests.get(url, cookies=cj)
        for chunk in read.iter_content(chunk_size=512):
            if chunk: fout.write(chunk)
        
        fout.close()
        
    fin.close()




def get_differential(condition, meta, data):
    """ Generate differential expression profile from meta data """
    
    # common conditions
    common = data.columns.intersection(meta.index).intersection(condition.index)
    
    N = common.shape[0]
    if N <= 1: return None
    
    #if N < data.shape[1] or N < meta.shape[0]: sys.stderr.write('Warning: not all samples are included\n')
    
    data = data.loc[:, common]
    meta = meta.loc[common]
    condition = condition.loc[common]
    
    meta_group = meta.groupby(condition)

    result_lst = []

    for condition, meta in meta_group:
        cnt_map = meta.value_counts()
        if len(cnt_map) <= 1 or 'Control' not in cnt_map.index: continue
        
        diff = data.loc[:, meta.index]
        diff = diff.groupby(meta.loc[diff.columns], axis=1).median()
        diff = diff.subtract(diff['Control'], axis=0).drop('Control', axis=1)
        
        rep_count = cnt_map.loc[diff.columns].apply(lambda v: min(v, cnt_map.loc['Control'])).astype(str)
        
        diff.columns += '@' + condition + ' rep ' + rep_count
        
        result_lst.append(diff)
    
    if len(result_lst) == 0:
        return None
    
    elif len(result_lst) == 1:
        # no need to merge
        return result_lst[0]
    
    else:
        # multiple conditions to one joint matrix on common genes
        return pandas.concat(result_lst, axis=1, join='inner')
    


def merge_sub_condition(result):
    """ Certain datasets may have sub conditions, merge them to generate a master differential files """
    info = pandas.DataFrame([v.split('@') for v in result.columns], index=result.columns, columns=['Treatment', 'Condition'])
    
    # .split(' join ')[0]
    info.loc[:, 'Count'] = info['Condition'].apply(lambda v: v.split(' rep ')[1]).astype(int)
    
    Condition = pandas.DataFrame(info['Condition'].apply(
        lambda v: pandas.Series(v.split(' rep ')[0].split('&'))
        ))
    
    Condition = Condition.drop(1, axis=1).apply(lambda v: '&'.join(v), axis=1)
    
    info.loc[:, 'Condition'] = info['Treatment'] + '@' + Condition
    
    cntmap = info['Count'].groupby(info['Condition']).sum()
    
    included = cntmap.index[cntmap > 1]
    if len(included) == 0: return None
    
    result.columns = info['Condition']
    result = result.groupby(result.columns, axis=1).median()
    
    return result.loc[:, included], cntmap



def simple_group(result):
    """ Ignore all conditions, just generate a merge file with conditions as biological replicates """
    
    result.columns = [v.split(' rep ')[0].split('@')[0].split('&')[0] for v in result.columns]
    
    cnt_map = result.columns.value_counts()
        
    if cnt_map.max() > 1:   
        result = result.groupby(result.columns, axis=1).median()
        result = result.loc[:, cnt_map.index[cnt_map > 1]]
            
        return result, cnt_map
    
    # still not cntmap with replicates merged
    return None



def process_curation_result_path(raw_path, output_path):
    """ Master function to analyze data annotations and generate differential expression profile """
    
    # get all meta annotations, each meta annotation will generate one differential profile
    meta_files = glob(os.path.join(raw_path, '*.meta.*'))
    
    for meta in meta_files:
        title = os.path.basename(meta)
        curator = title.split('.meta.')[1]
        
        #if title.split('.meta.')[0] not in ['GSE133968']: continue
        
        # find the matched expression matrices
        data_files = glob(re.sub('[.]meta.*', '.*.processed.gz', meta))
        if len(data_files) == 0: continue
        
        meta = pandas.read_csv(meta, sep='\t', index_col=0)
        
        # these two key columns must exist
        if 'Treatment' not in meta.columns or 'Condition' not in meta.columns:
            sys.stderr.write('Cannot find Treatment, Condition columns for %s\n' % title)
            continue
        
        meta = meta.loc[meta[['Condition', 'Treatment']].isnull().sum(axis=1) == 0]
        
        # two few samples annotated, no need to analyze such dataset
        if meta.shape[0] < 2:
            sys.stderr.write('Not sufficient samples for %s\n' % title)
            continue
        
        # at least some non-null values for each column
        meta = meta.loc[:, (~meta.isnull()).sum() > 0]    
        meta.fillna('', inplace=True)
        
        # does the current dataset have sub conditions?
        flag_sub_condition = ('Sub Condition' in meta.columns)
        
        # escape all separators: & @, because we will use these characters as field separators in column names
        meta = meta.astype(str).apply(lambda arr: arr.apply(lambda v: v.replace('&','-').replace('@','-')))
        
        # all possible first order combinations
        group_candidates = meta.columns.difference(['Condition', 'Sub Condition', 'Treatment'])
        
        groups = [[['Condition', 'Sub Condition', v], group_candidates.difference([v])] for v in group_candidates]
        
        groups.append([['Condition', 'Sub Condition'], group_candidates.difference(['Condition'])])
        
        for data_file in data_files:
            output = re.sub('[.]processed.gz', '.diff', data_file) + '.' + curator
            output = os.path.join(output_path, os.path.basename(output))
            
            data = pandas.read_csv(data_file, sep='\t', index_col=0)
            
            flag_write = False
            
            for condition, decorator in groups:
                print(condition, decorator)
                
                meta_sub = meta[[v for v in condition if v in meta.columns]]
                
                # make sure there is no empty space in every columns for group condition
                meta_sub = meta_sub.replace(r'^\s*$', numpy.nan, regex=True).dropna()
                
                # use escaped separator for split later
                condition = meta_sub.apply(lambda v: '&'.join(v.index + ':' + v), axis=1)
                
                treatment = copy.deepcopy(meta['Treatment'])
                
                if len(decorator) > 0:
                    # don't care the split of decorator
                    decorator = meta[decorator].apply(lambda v: '_'.join(v.index + ':' + v), axis=1)
                    
                    treatment += '&' + decorator  
                    treatment = treatment.apply(lambda v: v.rstrip('&').strip() if v.find('Control') < 0 else v.split('&')[0])
                
                result = get_differential(condition, treatment, data)
                
                # start at the first successful details
                if result is not None:
                    result = result.loc[(result == 0).mean(axis=1) < 1]
                    
                    if flag_sub_condition:
                        result.to_csv(output + '.sep.gz', sep='\t', index_label=False, compression='gzip')
                        
                        result = merge_sub_condition(result)
                        
                        if result is None:
                            sys.stderr.write('Failed merging sub conditions for %s\n' % data_file)
                            continue
                        else:
                            result, cntmap = result
                    
                    else:
                        # .split(' join ')[0]
                        info = pandas.DataFrame([v.split(' rep ') for v in result.columns], index=result.columns, columns = ['Condition', 'Count'])
                        
                        cntmap = info['Count'].astype(int)
                        
                        flag = (cntmap > 1).tolist()
                        
                        if sum(flag) > 0:
                            result = result.loc[:, flag]
                            result.columns = info.loc[flag, 'Condition']
                        else:
                            result.to_csv(output + '.sep.gz', sep='\t', index_label=False, compression='gzip')
                            result = simple_group(result)
                            
                            if result is None:
                                sys.stderr.write('Not sufficient replicates for %s.\n' % data_file)
                                continue
                            else:
                                result, cntmap = result
                    
                    result.to_csv(output, sep='\t', index_label=False)
                    cntmap.to_csv(output + '.cntmap', sep='\t', index_label=False, header=False)
                    
                    flag_write = True
                    
                    break
            
            if not flag_write: sys.stderr.write('None value for %s\n' % data_file)


def main():
    # download file list
    f = sys.argv[1]
    
    fpath = os.path.dirname(f)
    
    # download file path
    raw_path = os.path.join(fpath, 'raw')
    if not os.path.exists(raw_path): os.mkdir(raw_path)
    
    # differential output path
    diff_path = os.path.join(fpath, 'diff')
    if not os.path.exists(diff_path): os.mkdir(diff_path)
    
    # Step 1: download curated files from FDC server
    download_file_list(f, raw_path)
    
    # Step 2: generate differential expression files
    process_curation_result_path(raw_path, diff_path)
    
    return 0

if __name__ == '__main__': main()
