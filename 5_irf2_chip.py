#!/usr/bin/python


'''
The MIT License (MIT)

Copyright (c) 2018 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''


#Main method run script for processing of slam seq analysis from Muhar et al., 2018




#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys, os
# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = '/storage/cylin/bin/pipeline/'

sys.path.append(whereAmI)
sys.path.append(pipeline_dir)

import pipeline_dfci
import utils
import string
import numpy
import os
import re
from collections import defaultdict
import subprocess
#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = 'NIBR_YvsO_cyl'
genome ='hg19'
annotFile = '%s/annotation/%s_refseq.ucsc' % (pipeline_dir,genome)

#project folders
projectFolder = '/storage/cylin/grail/projects/%s' % (projectName) #PATH TO YOUR PROJECT FOLDER


projectFolder = utils.formatFolder(projectFolder,True)
#standard folder names
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)
metaRoseFolder = '%smeta_rose/' % (projectFolder)
roseFolder = '%srose/' % (projectFolder)
fastaFolder = '%sfasta/' % (projectFolder)
bedFolder = '%sbed/' % (projectFolder)
figuresFolder = '%sfigures/' % (projectFolder)
geneListFolder = '%sgeneListFolder/' % (projectFolder)
bedFolder = '%sbeds/' % (projectFolder)
signalFolder = '%ssignalTables/' % (projectFolder)
tableFolder = '%stables/' % (projectFolder)

#mask Files


#genomeDirectory #select your genome
#genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/'
#genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/hg19/Sequence/Chromosomes/'

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,roseFolder,fastaFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)


py27_path = '/storage/cylin/anaconda3/envs/py27_anaconda/bin/python'
#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================

#this project will utilize multiple datatables
#data tables are organized largely by type/system
#some data tables overlap for ease of analysis

#ChIP-Seq
chip_data_file = '%sdata_tables/NIBR_IRF2_CHIP_TABLE.txt' % (projectFolder)




#==========================================================================
#===========================MAIN METHOD====================================
#==========================================================================


def main():


    print('main analysis for project %s' % (projectName))

    print('changing directory to project folder')
    os.chdir(projectFolder)

    print('\n\n')
    print('#======================================================================')
    print('#======================I. LOADING DATA ANNOTATION======================')
    print('#======================================================================')
    print('\n\n')

    #This section sanity checks each data table and makes sure both bam and .bai files are accessible

    #for data file
    pipeline_dfci.summary(chip_data_file)


    print('\n\n')
    print('#======================================================================')
    print('#===========================II. CALLING MACS===========================')
    print('#======================================================================')
    print('\n\n')

    #pipeline_dfci.run_macs(chip_data_file,projectFolder,macsFolder,macsEnrichedFolder,wiggleFolder,useBackground=True)


    print('\n\n')
    print('#======================================================================')
    print('#=======================III. MERGING IRF2 REGIONS======================')
    print('#======================================================================')
    print('\n\n')

    #create a set of regions representing the intersect of peaks
    #filter out anything that overlaps a peak in the HA ctl


    def merge_regions():
        '''
        merges ha peaks to identify all overlapping peaks
        filters out anything overlapping the HA controls
        '''
        hk_dox_ha_1 = utils.importBoundRegion('%sHK_DOX_HA_1_peaks.bed' % (macsEnrichedFolder),'HK_DOX_HA_1')
        hk_dox_ha_2 = utils.importBoundRegion('%sHK_DOX_HA_2_peaks.bed' % (macsEnrichedFolder),'HK_DOX_HA_2')

        hk_dox_loci = hk_dox_ha_1.getLoci() + hk_dox_ha_2.getLoci()
        
        #control datasets
        hk_ctl_ha_1 = utils.importBoundRegion('%sHK_CTL_HA_1_peaks.bed' % (macsEnrichedFolder),'HK_CTL_HA_1')
        hk_ctl_ha_2 = utils.importBoundRegion('%sHK_CTL_HA_2_peaks.bed' % (macsEnrichedFolder),'HK_CTL_HA_2')

        hk_ctl_loci = hk_ctl_ha_1.getLoci() + hk_ctl_ha_2.getLoci()
        hk_ctl_lc = utils.LocusCollection(hk_ctl_loci)

        print(len(hk_dox_loci))
        stitched_lc = utils.LocusCollection(hk_dox_loci).stitchCollection()
        print(len(stitched_lc))
        filtered_loci = []
        for locus in stitched_lc.getLoci():
            if len(hk_dox_ha_1.getOverlap(locus)) > 0 and len(hk_dox_ha_2.getOverlap(locus)) > 0:
                if len(hk_ctl_lc.getOverlap(locus)) == 0:
                       filtered_loci.append(locus)
                

        print(len(filtered_loci))
        filtered_lc = utils.LocusCollection(filtered_loci)
        gff_path = '%sHG19_IRF2_HA_MERGED_FILTERED_CONSERVED_0_+0.gff' % (gffFolder)
        filtered_gff = utils.locusCollectionToGFF(filtered_lc)
        utils.unParseTable(filtered_gff,gff_path,'\t')
        
    #merge_regions()

    print('\n\n')
    print('#======================================================================')
    print('#======================IV. IDENTIFY ATAC OVERLAP REGIONS===============')
    print('#======================================================================')
    print('\n\n')

    # atac_bed_path = '%sHG19_combined_atac_-0_+0.bed' % (bedFolder)# all combined atac regions

    # atac_collection = utils.importBoundRegion(atac_bed_path,'HG19_combined_atac')
    # print(len(atac_collection))

    # #now filter the irf2 gff
    # irf2_gff_path = '%sHG19_IRF2_HA_MERGED_FILTERED_CONSERVED_0_+0.gff' % (gffFolder)
    
    # irf2_collection = utils.gffToLocusCollection(irf2_gff_path)
    # irf2_loci = irf2_collection.getLoci()
    
    # irf2_atac_loci = [locus for locus in irf2_loci if atac_collection.getOverlap(locus)]
    # print(len(irf2_atac_loci))
    # irf2_atac_collection=utils.LocusCollection(irf2_atac_loci)

    # irf2_atac_gff = utils.locusCollectionToGFF(irf2_atac_collection)
    # irf2_atac_gff_path = '%sHG19_IRF2_HA_MERGED_FILTERED_CONSERVED_ATAC_0_+0.gff' % (gffFolder)
    # utils.unParseTable(irf2_atac_gff,irf2_atac_gff_path,'\t')

    # # overlap with TSS
    # tss_gff_path = '%sHG19_TSS_ALL_-1000_+1000.gff' % (gffFolder)
    # tss_gff = utils.parseTable(tss_gff_path,'\t')
    # tss_collection = utils.gffToLocusCollection(tss_gff)
    
    # print('tss overlap w/ IRF2  atac peaks')
    # print(len([locus for locus in irf2_atac_loci if tss_collection.getOverlap(locus)]))
    # print(len(irf2_atac_loci))

    # #overlap w/ k27ac
    # k27ac_gff_path = '%sHG19_keratinocyte_combined_all_-0_+0.gff' % (gffFolder)
    # k27ac_gff = utils.parseTable(k27ac_gff_path,'\t')
    # k27ac_collection = utils.gffToLocusCollection(k27ac_gff)
    
    # print('k27ac overlap w/ IRF2  atac peaks')
    # print(len([locus for locus in irf2_atac_loci if k27ac_collection.getOverlap(locus)]))
    # print(len(irf2_atac_loci))



    print('\n\n')
    print('#======================================================================')
    print('#========================V. CALLING ROSE2 META=========================')
    print('#======================================================================')
    print('\n\n')

    def wrapRose2Meta(data_file,input_path,parent_folder,active_gene_path='',rank_list=[],control_list=[],analysis_name=''):
        '''
        quick wrapper for Rose2Meta
        '''
        dataDict = pipeline_dfci.loadDataTable(data_file)
        rank_string = ','.join([dataDict[name]['bam'] for name in rank_list])
        control_string = ','.join([dataDict[name]['bam'] for name in control_list])

        output_folder = utils.formatFolder('%s%s' % (parent_folder,analysis_name),True)    
        rose2_meta_cmd = '%s %sROSE2_META.py -g %s -i %s -r %s -c %s -n %s -o %s -s 0 -t 0 --mask %s' % (py27_path,pipeline_dir,genome,input_path,rank_string,control_string,analysis_name,output_folder,blacklist_path)

        all_enhancer_path = '%s%s_AllEnhancers.table.txt' % (output_folder,analysis_name)
        
        if active_gene_path != '':
            rose2_map_cmd = '%s %sROSE2_geneMapper.py -g %s -i %s -l %s' % (py27_path, pipeline_dir,genome,all_enhancer_path,active_gene_path)
        else:
            rose2_map_cmd = '%s %sROSE2_geneMapper.py -g %s -i %s' % (py27_path, pipeline_dir,genome,all_enhancer_path)


        rose_bash_path = '%s%s_rose2_meta.sh' % (parent_folder,analysis_name)
        rose_bash = open(rose_bash_path,'w')
        rose_bash.write('#!/usr/bin/python\n\n')
        rose_bash.write('#setting up bamliquidator\n')

        rose_bash.write('\n\n#ROSE2_CMD\n')
        rose_bash.write(rose2_meta_cmd +'\n')
        rose_bash.write(rose2_map_cmd +'\n')

        rose_bash.close()
        print('Wrote ROSE2 META CMD to %s' % (rose_bash_path))

    #use ROSE2 w/ -t 0 and -s 0 to quantify background subtracted AUC at all peaks


    # parent_folder = utils.formatFolder('%smeta_rose/' % (projectFolder),True)

    # blacklist_path = '/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Annotation/Masks/hg19_encode_blacklist.bed'
    # active_gene_path = '%sgeneListFolder/HG19_KERATINOCYTE_ACTIVE.txt' % (projectFolder)
    # #creating bam lists

    # rank_list = ['HK_DOX_HA_1','HK_DOX_HA_2']
    # control_list =['HK_DOX_WCE_1','HK_DOX_WCE_2']


    # #for all IRF2 HA
    # irf2_gff_path = '%sHG19_IRF2_HA_MERGED_FILTERED_CONSERVED_0_+0.gff' % (gffFolder)
    # analysis_name = 'IRF2_HA'
    # wrapRose2Meta(chip_data_file,irf2_gff_path,parent_folder,active_gene_path,rank_list,control_list,analysis_name)

    # irf2_atac_gff_path = '%sHG19_IRF2_HA_MERGED_FILTERED_CONSERVED_ATAC_0_+0.gff' % (gffFolder)
    # analysis_name = 'IRF2_HA_ATAC'
    # wrapRose2Meta(chip_data_file,irf2_atac_gff_path,parent_folder,active_gene_path,rank_list,control_list,analysis_name)

    
    print('\n\n')
    print('#======================================================================')
    print('#================VI. OVERLAPPING IRF2 W/ MOTIF PREDICTIONS=============')
    print('#======================================================================')
    print('\n\n')


    # #load up peaks
    # #irf2_atac_peaks
    # irf2_atac_gff_path = '%sHG19_IRF2_HA_MERGED_FILTERED_CONSERVED_ATAC_0_+0.gff' % (gffFolder)
    # irf2_atac_gff = utils.parseTable(irf2_atac_gff_path,'\t')
    # irf2_atac_loci = utils.gffToLocusCollection(irf2_atac_gff).getLoci()
    # irf2_atac_collection = utils.LocusCollection(irf2_atac_loci)
    # print(len(irf2_atac_loci))

    
    # irf2_edge_path = '%scrc_atac/keratinocyte_combined_all/keratinocyte_combined_all_EDGE_TABLE_signal_filtered_IRF2_EDGES.txt' % (projectFolder)

    # irf2_edge_table = utils.parseTable(irf2_edge_path,'\t')
    # print(len(irf2_edge_table))

    # irf2_confirmed_edges = []
    # irf2_edge_loci = []
    # for line in irf2_edge_table[1:]:
    #     chrom = line[1].split('(')[0]
    #     coords = [int(x) for x in line[1].split(':')[-1].split('-')]
    #     locus = utils.Locus(chrom,coords[0]-00,coords[1]+00,'.',line[0])
    #     if len(irf2_atac_collection.getOverlap(locus)) > 0:
    #         irf2_confirmed_edges.append(line)
    #     irf2_edge_loci.append(locus)
    # print(len(irf2_confirmed_edges))

    # irf2_confirmed_edge_path = '%scrc_atac/keratinocyte_combined_all/keratinocyte_combined_all_EDGE_TABLE_signal_filtered_IRF2_EDGES_CONFIRMED.txt' % (projectFolder)
    # utils.unParseTable(irf2_confirmed_edges,irf2_confirmed_edge_path,'\t')
    
    # irf2_edge_collection = utils.LocusCollection(irf2_edge_loci)
    # print(len(irf2_edge_collection))


    # overlap_count = 0
    # for locus in irf2_atac_loci:
    #     search_locus = utils.makeSearchLocus(locus,0,0)
    #     if len(irf2_edge_collection.getOverlap(search_locus)) >0:
    #         overlap_count+=1
    # print(overlap_count)

    print('\n\n')
    print('#======================================================================')
    print('#=================VII. RUNNING ENHANCER PROMOTER ON IRF2===============')
    print('#======================================================================')
    print('\n\n')


    def wrap_enhancer_promoter(dataFile,input_path,activity_path,analysis_name,names_list = [],useBackground=True):

        '''
        runs enhancer promoter on everybody with the conserved regions and union of active genes
        '''

        #hard coded paths
        tads_path ='%shESC_domains_hg19.bed' %(bedFolder)

        #setting the output folder
        ep_folder = utils.formatFolder('%senhancerPromoter/' % (projectFolder),True)

        dataDict = pipeline_dfci.loadDataTable(dataFile)
        if len(names_list) == 0:
            names_list = [name for name in dataDict.keys()]
            names_list.sort()

        bams_list = [dataDict[name]['bam'] for name in names_list]
        bams_string = ' '.join(bams_list)

        background_names = [dataDict[name]['background'] for name in names_list]
        background_list = [dataDict[background_name]['bam'] for background_name in background_names]
        background_string = ' '.join(background_list)


        ep_bash_path = '%s%s_enhancer_promoter.sh' % (ep_folder,analysis_name)
        ep_bash = open(ep_bash_path,'w')

        ep_bash.write('#!/usr/bin/bash\n\n\n')

        ep_bash.write('#enhancer promoter analysis for %s\n\n' % (analysis_name))

        if useBackground:
            python_cmd = 'python %senhancerPromoter.py -b %s -c %s -g %s -i %s -o %s -a %s --name %s --tads %s --top 2000\n\n' % (pipeline_dir,bams_string,background_string,genome.upper(),input_path,ep_folder,activity_path,analysis_name,tads_path)

            ep_bash.write(python_cmd)

        else:
            python_cmd = 'python %senhancerPromoter.py -b %s -g %s -i %s -o %s -a %s --name %s --tads %s --top 2000\n\n' % (pipeline_dir,bams_string,genome.upper(),input_path,ep_folder,activity_path,analysis_name,tads_path)

            ep_bash.write(python_cmd)

        ep_bash.close()

        return(ep_bash_path)


    # # blacklist_path = '/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Annotation/Masks/hg19_encode_blacklist.bed'
    # active_gene_path = '%sgeneListFolder/HG19_KERATINOCYTE_ACTIVE.txt' % (projectFolder)
    # irf2_atac_gff_path = '%sHG19_IRF2_HA_MERGED_FILTERED_CONSERVED_ATAC_0_+0.gff' % (gffFolder)
    # analysis_name = 'IRF2_HA_ATAC'
    # bam_list = ['HK_DOX_HA_1','HK_DOX_HA_2']    
    # wrap_enhancer_promoter(chip_data_file,irf2_atac_gff_path,active_gene_path,analysis_name,bam_list,useBackground=True)

    print('\n\n')
    print('#======================================================================')
    print('#==============VIII. FORMATTING THE HORRIFYING EXPRESSION TABLE========')
    print('#======================================================================')
    print('\n\n')

    # exp_path = '%sirf2_kd_rna_seq/single_counts_filtered_counts.txt' % (projectFolder)
    # sample_key_path = '%sirf2_kd_rna_seq/sample_key.txt' % (projectFolder)

    # sample_table = utils.parseTable(sample_key_path,'\t')
    # sample_list = [line[0] for line in sample_table[1:]]
    # print(sample_list)
    # exp_table = utils.parseTable(exp_path,'\t')

    # #for each gene make a dictionary
    # exp_dict = {}
    
    # #first fill out the dictionary by gene name
    # for line in exp_table[1:]:
    #     gene_name = line[3].replace('"','')
    #     exp_dict[gene_name] = {}

    # print(len(exp_dict.keys()))
    
    # for line in exp_table[1:]:
    #     gene_name = line[3].replace('"','')
    #     sample_ID = line[4].replace('"','')
    #     counts = line[2]
    #     exp_dict[gene_name][sample_ID] = counts

    # #make the formatted expression table
    # header = ['GENE_NAME'] + sample_list
    # exp_table_formatted = [header]
    # gene_list = exp_dict.keys()
    # gene_list.sort()
    # for gene in gene_list:
    #     exp_line = [gene] + [exp_dict[gene][sample_ID] for sample_ID in sample_list]
    #     exp_table_formatted.append(exp_line)

    # exp_table_formatted_path = '%sirf2_kd_rna_seq/irf2_expression_formatted.txt' % (projectFolder)
    # utils.unParseTable(exp_table_formatted,exp_table_formatted_path,'\t')
    
    # #with the exp dict we can make a nicer version of the gene table
    # gene_table_path = '%senhancerPromoter/IRF2_HA_ATAC/IRF2_HA_ATAC_GENE_TABLE.txt' % (projectFolder)
    # gene_table_formatted_path = '%senhancerPromoter/IRF2_HA_ATAC/IRF2_HA_ATAC_GENE_TABLE_FORMATTED.txt' % (projectFolder)

    # gene_table = utils.parseTable(gene_table_path,'\t')
    # gene_table_formatted = [gene_table[0] + ['IRF2_TOTAL_SIGNAL'] + header+ ['OLD_IRF2_KD_MEAN','OLD_CTL_MEAN','OLD_IRF2_VS_CTL','YOUNG_IRF2_KD_MEAN','YOUNG_CTL_MEAN','YOUNG_IRF2_VS_CTL']]
    # for line in gene_table[1:]:
    #     if float(line[1]) == 0.0 and float(line[2]) == 0.0:
    #         continue
    #     if exp_dict.has_key(line[0]) == False:
    #         continue
    #     gene = line[0]
    #     #where conditions are met
    #     old_kd_mean = numpy.mean([float(exp_dict[gene][x]) for x in ['IRF2_KD_OLD_1', 'IRF2_KD_OLD_2', 'IRF2_KD_OLD_3']])
    #     old_ctl_mean = numpy.mean([float(exp_dict[gene][x]) for x in ['CT_CRISPR_OLD_1', 'CT_CRISPR_OLD_2', 'CT_CRISPR_OLD_3']])
    #     old_fold = numpy.log2(old_kd_mean/old_ctl_mean)

    #     young_kd_mean = numpy.mean([float(exp_dict[gene][x]) for x in ['IRF2_KD_YOUNG_1', 'IRF2_KD_YOUNG_2', 'IRF2_KD_YOUNG_3']])
    #     young_ctl_mean = numpy.mean([float(exp_dict[gene][x]) for x in ['CT_CRISPR_YOUNG_1', 'CT_CRISPR_YOUNG_2', 'CT_CRISPR_YOUNG_3']])
    #     young_fold = numpy.log2(young_kd_mean/young_ctl_mean)

    #     exp_line = [gene] + [exp_dict[gene][sample_ID] for sample_ID in sample_list] + [round(x,4) for x in [old_kd_mean,old_ctl_mean,old_fold,young_kd_mean,young_ctl_mean,young_fold]]
    #     gene_table_formatted.append(line+[sum([float(x) for x in line[1:3]])] + exp_line)


    # utils.unParseTable(gene_table_formatted,gene_table_formatted_path,'\t')
        

    print('\n\n')
    print('#======================================================================')
    print('#=================IX. ANNOTATING IRF2 KD CLUSTERGRAM===================')
    print('#======================================================================')
    print('\n\n')

    #this little bit of python code is on the dropbox... need to move over

    print('\n\n')
    print('#======================================================================')
    print('#=======================X. PLOTTING FIGURE REGIONS ====================')
    print('#======================================================================')
    print('\n\n')


    figure_gff_path = '%sHG19_KERATINOCYTE_FIGURE_2_GENES.gff' % (gffFolder)
    plotName = 'IRF2_FIGURE_2_GENES'
    outputFolder = utils.formatFolder('%sgene_plot/IRF2/' % (projectFolder),True)
    pipeline_dfci.callBatchPlot(chip_data_file,figure_gff_path,plotName,outputFolder,namesList=['HK_DOX_HA_1','HK_DOX_HA_2'],uniform=True,bed ='',plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '',scaleFactorString ='')





#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
