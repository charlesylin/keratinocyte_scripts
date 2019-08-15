#!/usr/bin/python


'''
The MIT License (MIT)

Copyright (c) 2016 Charles Lin

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


#Main method run script for NIBR YvsO project

#1_map_enhancers.py runs enhancer calling and finds BRD4 subpeaks

#See README for additional information on downloading and installing dependencies

#Requires linlab pipeline set of utilities 
#

#Requires bamliquidator
#

#Requires 

#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys, os
# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = '/storage/cylin/home/cl6/pipeline/'

sys.path.append(whereAmI)
sys.path.append(pipeline_dir)

import pipeline_dfci
import utils
import string
import numpy
import os
import re
from collections import defaultdict
#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = 'NIBR_YvsO_cyl'
genome ='hg19'
annotFile = '%s/annotation/%s_refseq.ucsc' % (pipeline_dir,genome)

#project folders
projectFolder = '/storage/cylin/grail/projects/%s/' % (projectName) #PATH TO YOUR PROJECT FOLDER

#standard folder names
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)
metaRoseFolder = '%smeta_rose/' % (projectFolder)
fastaFolder = '%sfasta/' % (projectFolder)
bedFolder = '%sbed/' % (projectFolder)
figureCodeFolder = '%sfigureCode/' % (projectFolder)
figuresFolder = '%sfigures/' % (projectFolder)
geneListFolder = '%sgeneListFolder/' % (projectFolder)
bedFolder = '%sbeds/' % (projectFolder)
signalFolder = '%ssignalTables/' % (projectFolder)
tableFolder = '%stables/' % (projectFolder)
genePlotFolder = '%sgene_plot/' % (projectFolder)
#mask Files
maskFile ='%smasks/hg19_encode_blacklist.bed' % (projectFolder)

#genomeDirectory
genomeDirectory = '/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,fastaFolder,figureCodeFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder,genePlotFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)

py27_path = '/storage/cylin/anaconda3/envs/py27_anaconda/bin/python'
#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================

#this project will utilize multiple datatables
#data tables are organized largely by type/system
#some data tables overlap for ease of analysis

#ATAC-Seq
atac_dataFile = '%sdata_tables/NIBR_ATAC_TABLE_NEW_riesling.txt' % (projectFolder)

#ChIP-Seq
chip_dataFile = '%sdata_tables/NIBR_CHIP_TABLE.txt' % (projectFolder)

#RNA-Seq
rna_dataFile = '%sdata_tables/NIBR_RNA_TABLE.txt' % (projectFolder)


#==========================================================================
#===========================GLOBAL NAMES LISTS=============================
#==========================================================================

y_k27ac_list = ['Y_H3K27ac_1','Y_H3K27ac_2','Y_H3K27ac_3']
o_k27ac_list = ['O_H3K27ac_1','O_H3K27ac_2','O_H3K27ac_3']
k27ac_list = y_k27ac_list + o_k27ac_list

y_brd4_list = ['Y_BRD4_1','Y_BRD4_2']
o_brd4_list = ['O_BRD4_1','O_BRD4_2']

chip_list_no_input = [name for name in pipeline_dfci.loadDataTable(chip_dataFile).keys() if name.upper().count('INPUT') == 0]


#for some reason needed to re-run macs for that brd4
#pipeline_dfci.callMacsSlurm(chip_dataFile,macsFolder,namesList = ['O_BRD4_1'],overwrite=False,pvalue='1e-9',useBackground =True,launch=False)

#==========================================================================
#===========================MAIN METHOD====================================
#==========================================================================


def main():


    print('main analysis for MYCN project')

    print('changing directory to project folder')
    os.chdir(projectFolder)

    print('\n\n')
    print('#======================================================================')
    print('#========================I. CHECKING SEQ DATA==========================')
    print('#======================================================================')
    print('\n\n')

    #This section sanity checks each data table and makes sure both bam and .bai files are accessible

    #for ChIP-Seq
    #edit all of the data files to absolute path the

    #for ChIP-Seq
    pipeline_dfci.summary(chip_dataFile)

    #for RNA-Seq
    pipeline_dfci.summary(rna_dataFile)

    #for atac-seq
    pipeline_dfci.summary(atac_dataFile)

    print('\n\n')
    print('#======================================================================')
    print('#===============II. DEFINING H3K27AC ENHANCER LANDSCAPE================')
    print('#======================================================================')
    print('\n\n')

    # #for enhancers in young
    # analysis_name = 'young_h3k27ac'
    # y_enhancer_bashFileName,y_enhancer_region_map_path = define_enhancer_landscape(projectFolder,pipeline_dir,analysis_name,chip_dataFile,names_list= y_k27ac_list)

    # #runs only if no output detected
    # run_bash(y_enhancer_bashFileName,y_enhancer_region_map_path)


    # #for enhancers in old
    # analysis_name = 'old_h3k27ac'
    # o_enhancer_bashFileName,o_enhancer_region_map_path = define_enhancer_landscape(projectFolder,pipeline_dir,analysis_name,chip_dataFile,names_list= o_k27ac_list)

    # #runs only if no output detected
    # run_bash(o_enhancer_bashFileName,o_enhancer_region_map_path)


    print('\n\n')
    print('#======================================================================')
    print('#============IIa. DYNAMIC ANALYSIS OF YOUNG/OLD ENHANCERS==============')
    print('#======================================================================')
    print('\n\n')

    
    # #calls the dynamic enhancer code
    # analysis_name = 'young_old_h3k27ac_se'
    # name_1 = 'young_h3k27ac'
    # name_2 = 'old_h3k27ac'
    
    # rose_folder_1 = '%smeta_rose/%s/' % (projectFolder,name_1)
    # rose_folder_2 = '%smeta_rose/%s/' % (projectFolder,name_2)
    # rose_string = ','.join([rose_folder_1,rose_folder_2])

    # dynamic_folder = utils.formatFolder('%sdynamic_meta/' % (projectFolder),True)
    # output_folder = utils.formatFolder('%sdynamic_meta/%s/' % (projectFolder,analysis_name),True)
    
    # group_1_string = ','.join(y_k27ac_list)
    # group_2_string = ','.join(o_k27ac_list)

    # bash_path = '%s%s.sh' % (dynamic_folder,analysis_name)
    
    # bash_file = open(bash_path,'w')

    # bash_file.write('#!/usr/bin/bash\n\n')
    # bash_file.write('#SBATCH --mem=32000\n\n')


    # dynamic_cmd = '%s %sdynamicEnhancer_meta.py -g %s -d %s -r %s -o %s --group1 %s --group2 %s --name1 %s --name2 %s' % (py27_path, pipeline_dir, genome, chip_dataFile,rose_string,output_folder,group_1_string,group_2_string,name_1,name_2)
    
    # bash_file.write(dynamic_cmd)

    # bash_file.close()


    print('\n\n')
    print('#======================================================================')
    print('#===================III. RUNNING ENHANCER CLUSTERING===================')
    print('#======================================================================')
    print('\n\n')

    # analysis_name ='keratinocyte_se'
    # cluster_folder = utils.formatFolder('%sclustering' % (projectFolder),True)
    # cluster_rose_folder = utils.formatFolder('%sclustering_rose' % (projectFolder),True)
    # output_folder = utils.formatFolder('%s%s_clustering' % (cluster_folder,analysis_name),True)
    # names_string = ','.join(k27ac_list)
    # cluster_bash_path = '%s%s_clustering.sh' % (cluster_folder,analysis_name)
    # cluster_bash = open(cluster_bash_path,'w')
    # cluster_bash.write('#!/usr/bin/bash\n\n\n')
    # cluster_bash.write('#SBATCH --mem=16000\n\n\n')

    # cluster_cmd = '%s %sclusterEnhancer.py -d %s -i %s -r %s -o %s -e super -t 0 -n %s --mask %s' % (py27_path,pipeline_dir,chip_dataFile,names_string,cluster_rose_folder,output_folder,analysis_name,maskFile)
    # cluster_bash.write(cluster_cmd+'\n\n')
    # cluster_bash.close()

    #runs only if no output detected
    #run_bash(cluster_bash_path,output_path)

    # #for all enhancers
    # analysis_name ='keratinocyte_all'
    # cluster_folder = utils.formatFolder('%sclustering' % (projectFolder),True)
    # cluster_rose_folder = utils.formatFolder('%sclustering_rose' % (projectFolder),True)
    # output_folder = utils.formatFolder('%s%s_clustering' % (cluster_folder,analysis_name),True)
    # names_string = ','.join(k27ac_list)
    # cluster_bash_path = '%s%s_clustering.sh' % (cluster_folder,analysis_name)
    # cluster_bash = open(cluster_bash_path,'w')
    # cluster_bash.write('#!/usr/bin/bash\n\n\n')
    # cluster_bash.write('#SBATCH --mem=16000\n\n\n')

    # cluster_cmd = '%s %sclusterEnhancer.py -d %s -i %s -r %s -o %s -a -t 0 -n %s --mask %s' % (py27_path,pipeline_dir,chip_dataFile,names_string,cluster_rose_folder,output_folder,analysis_name,maskFile)
    # cluster_bash.write(cluster_cmd+'\n\n')
    # cluster_bash.close()

    # output_path = '%s%s_%s_clusterTable.txt' % (output_folder,genome.upper(),analysis_name)
    # print(output_path)
    # #runs only if no output detected
    # run_bash(cluster_bash_path,output_path)




    print('\n\n')
    print('#======================================================================')
    print('#================IV. MERGING KERATINOCYTE SE LANDSCAPES================')
    print('#======================================================================')
    print('\n\n')

    # #takes SE table output from previous section and quickly merges to make a gff
    # #also makes a gff individually

    # young_collection = utils.makeSECollection(y_enhancer_region_map_path,'young_h3k27ac')
    # old_collection = utils.makeSECollection(o_enhancer_region_map_path,'old_h3k27ac')

    # young_gff = utils.locusCollectionToGFF(young_collection)
    # old_gff = utils.locusCollectionToGFF(old_collection)

    # young_loci = young_collection.getLoci()
    # old_loci = old_collection.getLoci()

    # combined_loci = young_loci + old_loci 
    # combined_collection = utils.LocusCollection(combined_loci)
    # stitched_collection = combined_collection.stitchCollection()
    # combined_gff= utils.locusCollectionToGFF(stitched_collection)

    # #now write out some paths
    # young_gff_path = '%sHG19_young_h3k27ac_SE_-0_+0.gff' % (gffFolder)
    # old_gff_path = '%sHG19_old_h3k27ac_SE_-0_+0.gff' % (gffFolder)

    # #combined
    # combined_gff_path = '%sHG19_combined_h3k27ac_SE_-0_+0.gff' % (gffFolder)

    # #writing to disk
    # utils.unParseTable(young_gff,young_gff_path,'\t')
    # utils.unParseTable(old_gff,old_gff_path,'\t')
    # utils.unParseTable(combined_gff,combined_gff_path,'\t')


    print('\n\n')
    print('#======================================================================')
    print('#================V. MAKING BRD4 SUBPEAK BEDS USING MACS2===============')
    print('#======================================================================')
    print('\n\n')
    
    # #for now assumes pre-run macs2 will have to add that later
    # macs2Folder = '%smacs2Folder' % (projectFolder)
    # #for young
    # young_collection = mergeMacs2(macs2Folder,y_brd4_list)
    
    # #for old
    # old_collection = mergeMacs2(macs2Folder,o_brd4_list)


    # #making beds
    # young_bed = utils.locusCollectionToBed(young_collection)
    # old_bed = utils.locusCollectionToBed(old_collection)

    # #making combined collection
    # young_loci = young_collection.getLoci()
    # old_loci = old_collection.getLoci()

    # combined_loci = young_loci + old_loci 
    # combined_collection = utils.LocusCollection(combined_loci)
    # stitched_collection = combined_collection.stitchCollection()
    # combined_bed= utils.locusCollectionToBed(stitched_collection)

    # #now write out some paths
    # young_bed_path = '%sHG19_young_brd4_narrow_-0_+0.bed' % (bedFolder)
    # old_bed_path = '%sHG19_old_brd4_narrow_-0_+0.bed' % (bedFolder)

    # #combined
    # combined_bed_path = '%sHG19_combined_brd4_narrow_-0_+0.bed' % (bedFolder)

    # #writing to disk
    # utils.unParseTable(young_bed,young_bed_path,'\t')
    # utils.unParseTable(old_bed,old_bed_path,'\t')
    # utils.unParseTable(combined_bed,combined_bed_path,'\t')



    print('\n\n')
    print('#======================================================================')
    print('#===============VI. DEFINING ACTIVE GENES IN THE SYSTEM================')
    print('#======================================================================')
    print('\n\n')

    # # define active genes as FPKM > 10 in at least one sample
    # # and k27ac and brd4 present at the TSS

    # pipeline_dfci.makeGeneGFFs(annotFile,gffFolder,'HG19')

    # gffList = ['%sHG19_TSS_ALL_-1000_+1000.gff' % (gffFolder)]
    # cellTypeList = ['Y','O']
    
        
    # mapped_path = pipeline_dfci.mapEnrichedToGFF(chip_dataFile,'TSS',gffList,cellTypeList,macsEnrichedFolder,mappedEnrichedFolder,macs=True,namesList=chip_list_no_input,useBackground=True)
    # print(mapped_path)

    # mapped_path = '%sHG19_TSS_ALL_-1000_+1000/HG19_TSS_ALL_-1000_+1000_TSS.txt' % (mappedEnrichedFolder)
    # exp_path = '%s/cufflinks/NIBR_YvsO_cuffnorm/genes.fpkm_table' % (projectFolder)
    # activity_path = '%sHG19_KERATINOCYTE_ACTIVE.txt' % (geneListFolder)

    # makeActiveList(mapped_path,exp_path,annotFile,activity_path)
    

    print('\n\n')
    print('#======================================================================')
    print('#======================VII. RUNNING CIRCUITRY==========================')
    print('#======================================================================')
    print('\n\n')

    # #running circuitry on the consensus system
    # #creates a sbatch bash script

    # #for young
    # analysis_name = 'keratinocyte_young'
    # enhancer_path = '%smeta_rose/young_h3k27ac/young_h3k27ac_SuperEnhancers_ENHANCER_TO_GENE.txt' % (projectFolder)
    # subpeak_path = '%sbeds/HG19_young_brd4_narrow_-0_+0.bed' % (projectFolder)
    # activity_path = '%sgeneListFolder/HG19_KERATINOCYTE_ACTIVE.txt' % (projectFolder)

    # call_crc(analysis_name,enhancer_path,subpeak_path,activity_path)


    # #for old
    # analysis_name = 'keratinocyte_old'
    # enhancer_path = '%smeta_rose/old_h3k27ac/old_h3k27ac_SuperEnhancers_ENHANCER_TO_GENE.txt' % (projectFolder)
    # subpeak_path = '%sbeds/HG19_old_brd4_narrow_-0_+0.bed' % (projectFolder)
    # activity_path = '%sgeneListFolder/HG19_KERATINOCYTE_ACTIVE.txt' % (projectFolder)

    # call_crc(analysis_name,enhancer_path,subpeak_path,activity_path)

    # #for combined se
    # analysis_name = 'keratinocyte_combined'
    # enhancer_path = '%sclustering/keratinocyte_se_clustering/HG19_keratinocyte_se_clusterTable_ENHANCER_TO_GENE.txt' % (projectFolder)
    # subpeak_path = '%sbeds/HG19_combined_brd4_narrow_-0_+0.bed' % (projectFolder)
    # activity_path = '%sgeneListFolder/HG19_KERATINOCYTE_ACTIVE.txt' % (projectFolder)

    # call_crc(analysis_name,enhancer_path,subpeak_path,activity_path)


    # #for combined all enhancers
    # analysis_name = 'keratinocyte_combined_all'
    # enhancer_path = '%sclustering/keratinocyte_all_clustering/HG19_keratinocyte_all_clusterTable_ENHANCER_TO_GENE.txt' % (projectFolder)
    # subpeak_path = '%sbeds/HG19_combined_brd4_narrow_-0_+0.bed' % (projectFolder)
    # activity_path = '%sgeneListFolder/HG19_KERATINOCYTE_ACTIVE.txt' % (projectFolder)

    # call_crc(analysis_name,enhancer_path,subpeak_path,activity_path)



    print('\n\n')
    print('#======================================================================')
    print('#======================VIII. MAKE FORMATTED DEGREE TABLE===============')
    print('#======================================================================')
    print('\n\n')

    # young_degree_path = '%scrc/keratinocyte_young/keratinocyte_young_DEGREE_TABLE.txt' % (projectFolder)
    # young_degree_table  = utils.parseTable(young_degree_path,'\t')

    # old_degree_path = '%scrc/keratinocyte_old/keratinocyte_old_DEGREE_TABLE.txt' % (projectFolder)
    # old_degree_table  = utils.parseTable(old_degree_path,'\t')
    
    # tf_list = utils.uniquify([line[0] for line in young_degree_table[1:]] + [line[0] for line in old_degree_table[1:]])

    # degree_dict = {}
    # for tf in tf_list:

    #     degree_dict[tf] = [0.0,0.0,0.0,0.0]

    # young_count = len(young_degree_table)-1
    # for line in young_degree_table[1:]:
    #     tf = line[0]
    #     degree_dict[tf][0] = float(line[1])/young_count
    #     degree_dict[tf][1] = float(line[2])/young_count

    # old_count = len(old_degree_table)-1
    # for line in old_degree_table[1:]:
    #     tf = line[0]
    #     degree_dict[tf][2] = float(line[1])/old_count
    #     degree_dict[tf][3] = float(line[2])/old_count


    # tf_total_table = [['TF','YOUNG_IN','YOUNG_OUT','OLD_IN','OLD_OUT']]
    # tf_list.sort()
    # for tf in tf_list:
    #     tf_total_table.append([tf]+degree_dict[tf])

    # tf_total_path = '%scrc/keratinocyte_total_degree.txt' % (projectFolder)
    # utils.unParseTable(tf_total_table,tf_total_path,'\t')


    print('\n\n')
    print('#======================================================================')
    print('#======================VIII. PLOTTING TF REGIONS=======================')
    print('#======================================================================')
    print('\n\n')

    # #first make a gff of all enhancers
    # bed_string = '%scrc/keratinocyte_combined/keratinocyte_combined_all_subpeak.bed' % (projectFolder)
    # enhancer_tf_path = '%scrc/keratinocyte_combined/keratinocyte_combined_ENHANCER_TF_TABLE.txt' % (projectFolder)
    # enhancer_tf_table = utils.parseTable(enhancer_tf_path,'\t')
    # window = 50000
    # enhancer_tf_gff = []
    # for line in enhancer_tf_table[1:]:
    #     gff_line = [line[1],line[4],line[4],int(line[2])-window,int(line[3])+window,'','+','',line[4]]
    #     enhancer_tf_gff.append(gff_line)

    # figure_gff_path = '%sHG19_KERATINOCYTE_FIGURE_GENES.gff' % (gffFolder)
    # plot_prefix = 'HG19_KERATINOCYTE_FIGURE_1_GENES'
    # utils.unParseTable(enhancer_tf_gff,figure_gff_path,'\t')
    
    
    # #plot_keratinocyte_genes(plot_prefix,chip_dataFile,figure_gff_path,bed_string)    


    print('\n\n')
    print('#======================================================================')
    print('#====================IX. PLOTTING FIGURE REGIONS=======================')
    print('#======================================================================')
    print('\n\n')

    #plotting additional genes

    figure_2_gff = [['chr3','RBP1','RBP1',139232816,139260628,'','-','','RBP1'],
                    ['chr4','SMAD1','SMAD1',146394843,146483973,'','+','','SMAD1'],
                    ['chr4','UBA6','UBA6',68526969,68600569,'','-','','UBA6'],
                    ['chr15','CLN6','CLN6',68508894,68523908,'','-','','CLN6'],
                    ['chr3','XPC','XPC',14181835,14222159,'','-','','XPC'],
                    ['chr11','DCPS','DCPS',126170837,126196119,'','+','','DCPS'],
                    ['chr5','IRX2','IRX2',2735144,2765154,'','-','','IRX2'],
                    ['chr8','SNAI2_WIDE','SNAI2_WIDE',49746215,49846374,'','-','','SNAI2_WIDE'],
                    ['chr8','SNAI2_CLOSE','SNAI2_CLOSE',49826529,49841529,'','-','','SNAI2_CLOSE'],
                    ['chr12','MMP17','MMP17',132308545,132323642,'','+','','MMP17'],
                    ['chr15','VPS18','VPS18',41182805,41212811,'','+','','VPS18'],
                    ['chr4','ZFP42','ZFP42',188914912,188946979,'','+','','ZFP42'],
                    ['chrX','SAT1','SAT1',23794388,23827662,'','+','','SAT1'],
                    ['chr6','MAPK13','MAPK13',36096634,36101634,'','+','MAPK13'],
                    ['chr6','PSMB9','PSMB9',32800513,32831638,'','+','','PSMB9'],
                    ['chr6','TAP1','TAP1',32816769,32826769,'','+','','TAP1'],
                    ['chr11','TRIM44','TRIM44',35676627,35849111,'','+','','TRIM44'],
                    ['chr15','B2M','B2M',44996447,45012562,'','+','','B2M'],

                ]
    figure_gff_path = '%sHG19_KERATINOCYTE_FIGURE_2_GENES.gff' % (gffFolder)
    utils.unParseTable(figure_2_gff,figure_gff_path,'\t')


    bed_list = ['%scrc_atac/keratinocyte_combined_all/motif_beds/IRF2_motifs.bed' % (projectFolder),
                '%scrc_atac/keratinocyte_combined_all/motif_beds/SNAI2_motifs.bed' % (projectFolder),
                '%scrc_atac/keratinocyte_combined_all/keratinocyte_combined_all_all_subpeak.bed' % (projectFolder),
                ]
    
    bed_string = ','.join(bed_list)
    plot_prefix = 'HG19_KERATINOCYTE_FIGURE_2_GENES'
    plot_keratinocyte_genes(plot_prefix,chip_dataFile,figure_gff_path,bed_string)    


    #plotting additional genes atac

    bed_list = ['%scrc_atac/keratinocyte_combined_all/motif_beds/IRF2_motifs.bed' % (projectFolder),
                '%scrc_atac/keratinocyte_combined_all/motif_beds/SNAI2_motifs.bed' % (projectFolder),
                '%scrc_atac/keratinocyte_combined_all/keratinocyte_combined_all_all_subpeak.bed' % (projectFolder),
                ]
    
    bed_string = ','.join(bed_list)

    atac_dataFile_riesling = '%sdata_tables/NIBR_ATAC_TABLE_NEW_riesling.txt' % (projectFolder)
    plot_prefix = 'HG19_KERATINOCYTE_FIGURE_2_GENES_ATAC'
    plot_keratinocyte_atac(plot_prefix,atac_dataFile_riesling,figure_gff_path,bed_string)  
    





    print('\n\n')
    print('#======================================================================')
    print('#==================X. RUNNING ENHANCER PROMOTER QUANT FOR ALL==========')
    print('#======================================================================')
    print('\n\n')

    #quanitfies BRD4 at individual genes in Y vs O

    # first needs a gff of all of the enhancers
    # all_enhancer_path = '%scrc/keratinocyte_combined_all/keratinocyte_combined_all_ENHANCER_TABLE.txt' % (projectFolder)

    # all_enhancer_gff = []
    # all_enhancer_gff_path = '%sHG19_keratinocyte_combined_all_-0_+0.gff' % (gffFolder)

    # all_enhancer_table = utils.parseTable(all_enhancer_path,'\t')    
    # for line in all_enhancer_table[1:]:
    #     gff_line = [line[1],line[0],line[0],line[2],line[3],'','.','',line[0]]
    #     all_enhancer_gff.append(gff_line)

    # utils.unParseTable(all_enhancer_gff,all_enhancer_gff_path,'\t')


    # #for young
    # input_path = '%sHG19_keratinocyte_combined_all_-0_+0.gff' % (gffFolder)
    # activity_path = '%sgeneListFolder/HG19_KERATINOCYTE_ACTIVE.txt' % (projectFolder)
    # analysis_name = 'keratinocyte_combined_all_young'
    # y_brd4_list = ['Y_BRD4_1','Y_BRD4_2']
    # output_folder = utils.formatFolder('%senhancerPromoter/%s' % (projectFolder,analysis_name),True)
    # enhancer_promoter_bash,enhancer_promoter_output = wrap_enhancer_promoter(chip_dataFile,input_path,activity_path,analysis_name,y_brd4_list)
    # print(enhancer_promoter_output)
    # run_bash(enhancer_promoter_bash,enhancer_promoter_output)

    # #for young_h3k27ac
    # input_path = '%sHG19_keratinocyte_combined_all_-0_+0.gff' % (gffFolder)
    # activity_path = '%sgeneListFolder/HG19_KERATINOCYTE_ACTIVE.txt' % (projectFolder)
    # analysis_name = 'keratinocyte_combined_all_young_h3k27ac'
    # y_h3k27ac_list = ['Y_H3K27ac_1','Y_H3K27ac_2','Y_H3K27ac_3']
    # output_folder = utils.formatFolder('%senhancerPromoter/%s' % (projectFolder,analysis_name),True)
    # enhancer_promoter_bash,enhancer_promoter_output = wrap_enhancer_promoter(chip_dataFile,input_path,activity_path,analysis_name,y_h3k27ac_list)
    # print(enhancer_promoter_output)
    # run_bash(enhancer_promoter_bash,enhancer_promoter_output)


    # #for old
    # input_path = '%sHG19_keratinocyte_combined_all_-0_+0.gff' % (gffFolder)
    # activity_path = '%sgeneListFolder/HG19_KERATINOCYTE_ACTIVE.txt' % (projectFolder)
    # analysis_name = 'keratinocyte_combined_all_old'
    # o_brd4_list = ['O_BRD4_1','O_BRD4_2']
    # output_folder = utils.formatFolder('%senhancerPromoter/%s' % (projectFolder,analysis_name),True)
    # enhancer_promoter_bash,enhancer_promoter_output = wrap_enhancer_promoter(chip_dataFile,input_path,activity_path,analysis_name,o_brd4_list)
    # run_bash(enhancer_promoter_bash,enhancer_promoter_output)


    # #for old h3k27ac
    # input_path = '%sHG19_keratinocyte_combined_all_-0_+0.gff' % (gffFolder)
    # activity_path = '%sgeneListFolder/HG19_KERATINOCYTE_ACTIVE.txt' % (projectFolder)
    # analysis_name = 'keratinocyte_combined_all_old_h3k27ac'
    # o_h3k27ac_list = ['O_H3K27ac_1','O_H3K27ac_2','O_H3K27ac_3']
    # output_folder = utils.formatFolder('%senhancerPromoter/%s' % (projectFolder,analysis_name),True)
    # enhancer_promoter_bash,enhancer_promoter_output = wrap_enhancer_promoter(chip_dataFile,input_path,activity_path,analysis_name,o_h3k27ac_list)
    # run_bash(enhancer_promoter_bash,enhancer_promoter_output)



    print('\n\n')
    print('#======================================================================')
    print('#===============XI. DIFFERENTIAL ENHANCER PROMOTER ANALYSIS============')
    print('#======================================================================')
    print('\n\n')

    # #create a delta waterfall of the enhancer promoter output between young and old
    # #takes a pair of enhancer promoter gene tables and calls a simple R script to make a waterfall

    # #correlates to changes in expression
    # exp_path = '%scufflinks/NIBR_YvsO_cuffnorm/genes.fpkm_table' % (projectFolder)

    # name_1 = 'keratinocyte_combined_all_old'
    # ep_gene_path_1 = '%senhancerPromoter/%s/%s_GENE_TABLE.txt' % (projectFolder,name_1,name_1)

    # name_2 = 'keratinocyte_combined_all_young'
    # ep_gene_path_2 = '%senhancerPromoter/%s/%s_GENE_TABLE.txt' % (projectFolder,name_2,name_2)

    # h3k27ac_name_1 = 'keratinocyte_combined_all_old_h3k27ac'
    # h3k27ac_ep_gene_path_1 = '%senhancerPromoter/%s/%s_GENE_TABLE.txt' % (projectFolder,h3k27ac_name_1,h3k27ac_name_1)

    # h3k27ac_name_2 = 'keratinocyte_combined_all_young_h3k27ac'
    # h3k27ac_ep_gene_path_2 = '%senhancerPromoter/%s/%s_GENE_TABLE.txt' % (projectFolder,h3k27ac_name_2,h3k27ac_name_2)

    # output_path = '%sfigures/%s_%s_enhancer_promoter_gene_delta.pdf' % (projectFolder,name_1,name_2)

    # waterfall_script_path = '%sr_scripts/enhancer_promoter_waterfall.R' % (projectFolder)
    
    # waterfall_bash_script = '%sr_scripts/enhancer_promoter_waterfall_run.sh' % (projectFolder)
    # waterfall_bash = open(waterfall_bash_script,'w')
    
    # waterfall_bash.write('#!/usr/bin/bash\n\n')
    
    
    # r_cmd = 'Rscript %s %s %s %s %s %s %s %s %s' % (waterfall_script_path,name_1,name_2,ep_gene_path_1,ep_gene_path_2,h3k27ac_ep_gene_path_1,h3k27ac_ep_gene_path_2,exp_path,output_path)

    # waterfall_bash.write(r_cmd)
    # waterfall_bash.close()
    # print(r_cmd)
    # #os.system(r_cmd)





    print('\n\n')
    print('#======================================================================')
    print('#==================XII. DIFFERENTIAL OUT DEGREE ANALYSIS===============')
    print('#======================================================================')
    print('\n\n')

    #this is deprecated see 4_run_crc_atac.py

    #two ways to look for circuitry changes between Y and O
    #1 anchor on SE associated TFs, look at all of their assigned subpeaks
    #calcualte the delta brd4 for each
    #first need to map regions to the brd4 subpeaks
    #then build dictionaries


    # crc_folder = '%scrc/keratinocyte_combined/' % (projectFolder)
    # analysis_name = 'keratinocyte_combined'
    # tf_motif_brd4_delta(crc_folder,chip_dataFile,analysis_name,y_brd4_list,o_brd4_list)


    # crc_folder = '%scrc/keratinocyte_combined_all/' % (projectFolder)
    # analysis_name = 'keratinocyte_combined_all'
    # tf_motif_brd4_delta(crc_folder,chip_dataFile,analysis_name,y_brd4_list,o_brd4_list)



    # # #by groups of tfs
    # analysis_name = 'keratinocyte_combined'
    # crc_folder = '%scrc/%s/' % (projectFolder,analysis_name)
    


    # # tf_list = ['ELK3','FOSL1','VDR']
    # # tf_group_name = 'UP'
    # # all_3_out = tf_group_brd4_delta(crc_folder,chip_dataFile,analysis_name,tf_list,tf_group_name,y_brd4_list,o_brd4_list)


    # tf_list = []
    # tf_group_name = 'all'
    # all_out = tf_group_brd4_delta(crc_folder,chip_dataFile,analysis_name,tf_list,tf_group_name,y_brd4_list,o_brd4_list)




    # # #by groups of tfs
    # analysis_name = 'keratinocyte_combined_all'
    # crc_folder = '%scrc/%s/' % (projectFolder,analysis_name)
    



    # tf_list = []
    # tf_group_name = 'all'
    # all_out = tf_group_brd4_delta(crc_folder,chip_dataFile,analysis_name,tf_list,tf_group_name,y_brd4_list,o_brd4_list)








#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


#specific functions written for this analysis


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~RUN AND CHECK UTILITY~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#simple helper to run and check a bash script

def run_bash(bash_path,output_path,maxWait=30):
    '''
    runs a bash script and waits up to N minutes
    '''

    if not utils.checkOutput(output_path,0,0):

        print('running bash script %s' % (bash_path))
        os.system('bash %s' % (bash_path))
        if utils.checkOutput(output_path,1,30):
            print('run completed, output detected for %s at %s' % (bash_path,output_path))
    else:
        print('found prior output for %s at %s' % (bash_path,output_path))

            

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~DEFINING KERATINOCYTE H3K27AC ENHANCER LANDSCAPE~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def define_enhancer_landscape(projectFolder,pipeline_dir,analysis_name,chip_dataFile,names_list= []):

    '''
    defines the enhancer baseline using H3K27ac chips 
    enhancers defined using auto optimized stitching of nearby regions
    w/ a 2.5kb tss exclusion
    uses the meta rose code and writes out a .sh file for reproducibility
    '''

    #For H3K27AC
    #with TSS exclusion and auto stitching

    dataDict = pipeline_dfci.loadDataTable(chip_dataFile)

    if len(names_list) == 0:
        names_list = [name for name in dataDict.keys() if name.count('H3K27AC') == 1]

    bamFileList = [dataDict[name]['bam'] for name in names_list]
    bamString = string.join(bamFileList,',')

    controlBams = [dataDict[name]['background'] for name in names_list]
    controlFileList = [dataDict[name]['bam'] for name in controlBams]
    controlBamString = string.join(controlFileList,',')

    bedFileList = [macsEnrichedFolder + dataDict[name]['enrichedMacs'] for name in names_list]
    bedString = string.join(bedFileList,',')

    roseFolder = '%smeta_rose/' % (projectFolder)
    roseFolder = utils.formatFolder(roseFolder,True)

    outputFolder = '%s%s/' % (roseFolder,analysis_name)
    bashFileName = '%s%s_meta_rose.sh' % (roseFolder,analysis_name)

    bashFile = open(bashFileName,'w')
    bashFile.write('#!/usr/bin/bash\n\n')
    bashFile.write('#SBATCH --mem=16000\n\n') #set memory for slurm
    bashFile.write('cd %s\n' % (pipeline_dir))

    metaRoseCmd = 'python %sROSE2_META.py -g hg19 -i %s -r %s -c %s -o %s -n %s -t 0 --mask %s' % (pipeline_dir,bedString,bamString,controlBamString,outputFolder,analysis_name,maskFile)

    bashFile.write(metaRoseCmd + '\n')
    bashFile.close()


    #the 4KB parameter is 
    region_map_path = '%s%s/%s_SuperEnhancers.table.txt' % (roseFolder,analysis_name,analysis_name)
    return bashFileName,region_map_path


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~MAKING SUBPEAKS FROM MACS2 FILES~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def mergeMacs2(macs2Folder,names_list):

    '''
    merges a bunch of peaks across dataset names
    spits out as a locus collection
    '''
    macs2Folder = utils.formatFolder(macs2Folder,False)

    bed_loci = []
    for name in names_list:
        macs2_peak_path = '%s%s/%s_peaks.narrowPeak' % (macs2Folder,name,name)
        macs2_peaks = utils.parseTable(macs2_peak_path,'\t')
        peaks_loci = []
        for line in macs2_peaks:
            peaks_loci.append(utils.Locus(line[0],line[1],line[2],'.',line[3]))

        bed_loci += peaks_loci

    bed_collection = utils.LocusCollection(bed_loci)
    stitched_collection = bed_collection.stitchCollection()

    return stitched_collection
    

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~DEFINING ACTIVE GENE LISTS~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def makeActiveList(mapped_path,exp_path,annotFile,output=''):
    
    '''
    makes a list of active genes
    10fpkm expression cutoff
    bound by k27ac and brd4 in at least 1 dataset
    '''


    exp_table =utils.parseTable(exp_path,'\t') #BC1 = young BC=old
    
    exp_dict= defaultdict(float)
    for line in exp_table[1:]:
        exp_vector = [float(x) for x in line[1:]]
        exp_dict[line[0]] = max(exp_vector)
    
    startDict = utils.makeStartDict(annotFile)

    #now make the active gene list first by chromatin
    mapped_table = utils.parseTable(mapped_path,'\t')
    active_table = []
    for i in range(1,len(mapped_table)):
        line = mapped_table[i]
        #check old
        
        if max([int(line[i]) for i in [2,3]]) and max([int(line[i]) for i in [4,5,6]]) == 1:
            old_active = True
        else:
            old_active = False

        if max([int(line[i]) for i in [7,8]]) and max([int(line[i]) for i in [9,10,11]]) == 1:
            young_active = True
        else:
            young_active = False
    
        if old_active or young_active:
            chromatin_active = True
        
        gene_name = startDict[line[1]]['name']
        
        if exp_dict[gene_name] > 10.0:
            active_table.append([line[1],gene_name])

    if len(output) == 0:
        return active_table
    else:
        utils.unParseTable(active_table,output,'\t')
        return output


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~CALLING CRC ANALYSIS~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def call_crc(analysis_name,enhancer_path,subpeak_path,activity_path):

    '''
    runs crc
    '''


    crc_folder = utils.formatFolder('%scrc' % (projectFolder),True)

    output_folder = utils.formatFolder('%s%s' % (crc_folder,analysis_name),True)

    crc_cmd = '%s %scrc/CRC3.py -e %s -g HG19 -o %s -n %s -s %s -a %s' % (py27_path,pipeline_dir,enhancer_path,output_folder,analysis_name,subpeak_path,activity_path)
    crc_bash_path = '%s%s_crc.sh' % (crc_folder,analysis_name)

    crc_bash = open(crc_bash_path,'w')
    crc_bash.write('#!/usr/bin/bash\n\n')
    crc_bash.write('#SBATCH --mem=16000\n\n')
    crc_bash.write(crc_cmd +'\n')
    crc_bash.close()

    print('wrote crc command for %s to %s' % (analysis_name,crc_bash_path))
    return crc_bash_path
        

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~PLOTTING METAS ACROSS NB~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def plot_keratinocyte_genes(plot_prefix,chip_dataFile,figure_gff_path,bed_string):

    '''
    plots all varieties and iterations of tracks for shep on data
    '''


    #first establish the plot folder
    plotFolder = utils.formatFolder('%sKERATINOCYTE_ALL/' % (genePlotFolder),True)
    #plot_prefix = 'HG19_KERATINOCYTE_GENES'
    
    #we also have to set the extension properly between datasets

    #go by data file
    dataDict = pipeline_dfci.loadDataTable(chip_dataFile)
    names_list = dataDict.keys()
    #initial check for consistency of read lengths
    # for name in names_list:
    #     bam = utils.Bam(dataDict[name]['bam'])
    #     read_length = bam.getReadLengths()[0]
    #     bam_extension = 200-read_length
    #     print('For dataset %s in %s using an extension of %s' % (name,chip_dataFile,bam_extension))

    bam = utils.Bam(dataDict[names_list[0]]['bam'])
    read_length = bam.getReadLengths()[0]
    bam_extension = 200-read_length
    print('For datasets in %s using an extension of %s' % (chip_dataFile,bam_extension))
    
    #first do individuals
    for plot_group in ['BRD4','H3K27AC']:
        plotList = [name for name in dataDict.keys() if name.upper().count(plot_group) > 0 and name.upper().count('INPUT') == 0]
        plotName = '%s_INDIVIDUAL_%s' % (plot_prefix,plot_group)
        print(plotName)
        pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '')

    #now as metas
    plotList = y_k27ac_list + o_k27ac_list 
    groupList = ['YOUNG_H3K27AC']*3 + ['OLD_H3K27AC']*3 
    groupString = ','.join(groupList)

    plotName = '%s_H3K27AC_META_RELATIVE' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')


    plotName = '%s_H3K27AC_META_UNIFORM' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')


    #now as metas for brd4
    plotList = y_brd4_list + o_brd4_list
    groupList = ['YOUNG_BRD4']*2 + ['OLD_BRD4']*2
    groupString = ','.join(groupList)

    plotName = '%s_BRD4_META_RELATIVE' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')


    plotName = '%s_BRD4_META_UNIFORM' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')

    #now metas for all
    plotList =  y_k27ac_list + o_k27ac_list +y_brd4_list + o_brd4_list
    groupList = ['YOUNG_H3K27AC']*3 + ['OLD_H3K27AC']*3  +['YOUNG_BRD4']*2 + ['OLD_BRD4']*2
    groupString = ','.join(groupList)

    plotName = '%s_ALL_META_RELATIVE' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')


    plotName = '%s_ALL_META_UNIFORM' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')





def plot_keratinocyte_atac(plot_prefix,atac_dataFile,figure_gff_path,bed_string):

    '''
    plots all varieties and iterations of tracks for shep on data
    '''


    #first establish the plot folder
    plotFolder = utils.formatFolder('%sKERATINOCYTE_ATAC/' % (genePlotFolder),True)
    #plot_prefix = 'HG19_KERATINOCYTE_GENES'
    
    #we also have to set the extension properly between datasets

    #go by data file
    dataDict = pipeline_dfci.loadDataTable(atac_dataFile)
    names_list = dataDict.keys()
    #initial check for consistency of read lengths
    # for name in names_list:
    #     bam = utils.Bam(dataDict[name]['bam'])
    #     read_length = bam.getReadLengths()[0]
    #     bam_extension = 200-read_length
    #     print('For dataset %s in %s using an extension of %s' % (name,atac_dataFile,bam_extension))

    bam = utils.Bam(dataDict[names_list[0]]['bam'])
    read_length = bam.getReadLengths()[0]
    bam_extension = 200-read_length
    print('For datasets in %s using an extension of %s' % (atac_dataFile,bam_extension))
    
    #first do individuals
    for plot_group in ['Y','O']:
        plotList = [name for name in dataDict.keys() if name.upper().count(plot_group) > 0 and name.upper().count('INPUT') == 0]
        plotName = '%s_INDIVIDUAL_%s' % (plot_prefix,plot_group)
        print(plotName)
        pipeline_dfci.callBatchPlot(atac_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '')

    y_list = [name for name in dataDict.keys() if name.upper().count('Y') >0]
    o_list = [name for name in dataDict.keys() if name.upper().count('O') >0]
    #now as metas
    plotList = y_list + o_list 
    groupList = ['YOUNG']*len(y_list) + ['OLD']*len(o_list)
    groupString = ','.join(groupList)

    plotName = '%s_META_RELATIVE' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(atac_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')


    plotName = '%s_META_UNIFORM' % (plot_prefix)    
    pipeline_dfci.callBatchPlot(atac_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~RUNNING ENHANCER PROMOTER ANALYSIS~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def wrap_enhancer_promoter(dataFile,input_path,activity_path,analysis_name,names_list,useBackground=True,top=5000):

    '''
    runs enhancer promoter on everybody with the conserved regions and union of active genes
    '''
    
    #hard coded paths
    tads_path ='%shESC_domains_hg19.bed' %(bedFolder)

    #setting the output folder
    ep_folder = utils.formatFolder('%senhancerPromoter/' % (projectFolder),True)
    


    dataDict = pipeline_dfci.loadDataTable(dataFile)

    bams_list = [dataDict[name]['bam'] for name in names_list]
    bams_string = ' '.join(bams_list)
    
    background_names = [dataDict[name]['background'] for name in names_list]
    background_list = [dataDict[background_name]['bam'] for background_name in background_names]
    background_string = ' '.join(background_list)


    ep_bash_path = '%s%s_enhancer_promoter.sh' % (ep_folder,analysis_name)
    ep_bash = open(ep_bash_path,'w')

    ep_bash.write('#!/usr/bin/bash\n\n\n')
    ep_bash.write('#SBATCH --mem=16000\n\n')
    
    ep_bash.write('#enhancer promoter analysis for %s\n\n' % (analysis_name))

    ep_bash.write('cd %s%s/\n\n' % (ep_folder,analysis_name))
    if useBackground:

        python_cmd = 'python %senhancerPromoter.py -b %s -c %s -g %s -i %s -o %s -a %s --name %s --tads %s --top %s\n\n' % (pipeline_dir,bams_string,background_string,genome.upper(),input_path,ep_folder,activity_path,analysis_name,tads_path,top)

        ep_bash.write(python_cmd)
    else:

        python_cmd = 'python %senhancerPromoter.py -b %s -g %s -i %s -o %s -a %s --name %s --tads %s --top %s\n\n' % (pipeline_dir,bams_string,genome.upper(),input_path,ep_folder,activity_path,analysis_name,tads_path,top)

        ep_bash.write(python_cmd)

    ep_bash.close()
    
    #generate an anticipated output to check for completeness
    output_path = '%s%s/%s_TOP_%s_ORDERED.txt' % (ep_folder,analysis_name,analysis_name,top)

    return(ep_bash_path,output_path)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~CALCULATING CHANGES IN BRD4 OUT DEGREE BY TF~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def tf_motif_brd4_delta(crc_folder,chip_dataFile,analysis_name,y_brd4_list,o_brd4_list):

    '''
    calculates changes in brd4 out degree at each predicted motif occurrence
    '''
    crc_folder = utils.formatFolder(crc_folder,False)
    subpeak_bed_path = '%s%s_all_subpeak.bed' % (crc_folder,analysis_name)
    #direct the output to the crc folder
    signal_path = '%s%s_all_subpeak_signal.txt' % (crc_folder,analysis_name)
    
    #set final output
    output_path = '%s%s_brd4_subpeak_delta_out.txt' % (crc_folder,analysis_name)

    
    all_brd4_list = y_brd4_list + o_brd4_list
    pipeline_dfci.map_regions(chip_dataFile,[subpeak_bed_path],mappedFolder,signalFolder,all_brd4_list,True,signal_path)

    #now bring in the signal table as a dictionary using the locus line as the id
    print('making log2 y vs o signal dict at subpeaks')
    signal_table = utils.parseTable(signal_path,'\t')
    signal_dict = defaultdict(float)
    
    #figure out columns for young and old
    o_columns = [signal_table[0].index(name) for name in o_brd4_list]
    y_columns = [signal_table[0].index(name) for name in y_brd4_list]
    for line in signal_table[1:]:
        o_signal = numpy.mean([float(line[col]) for col in o_columns])
        y_signal = numpy.mean([float(line[col]) for col in y_columns])
        delta = numpy.log2(y_signal/o_signal)
        signal_dict[line[1]] = delta

    #loading subpeak collection
    subpeak_table = utils.parseTable(subpeak_bed_path,'\t')
    subpeak_loci = []

    for line in subpeak_table:
        line_locus = utils.Locus(line[0],int(line[1]),int(line[2]),'.')
        subpeak_loci.append(line_locus)

        
    #next create a dictionary of all tfs in the system
    tf_table_path = '%s%s_GENE_TF_TABLE.txt' % (crc_folder,analysis_name)
    
    tf_table = utils.parseTable(tf_table_path,'\t')
    tf_list = utils.uniquify([line[0] for line in tf_table[1:]])
    tf_list.sort()
    
    print('finding motifs for tfs')
    print(tf_list)


    out_degree_delta_dict = defaultdict(list)
    for tf_name in tf_list:
        print(tf_name)
        motif_bed_path = '%smotif_beds/%s_motifs.bed' % (crc_folder,tf_name)
        motif_bed = utils.parseTable(motif_bed_path,'\t')
        motif_loci = []
        for line in motif_bed[1:]:
            motif_locus = utils.Locus(line[0],line[1],line[2],'.')
            motif_loci.append(motif_locus)

        motif_collection = utils.LocusCollection(motif_loci)

        #now iterate through the subpeak loci
        for subpeak_locus in subpeak_loci:
            if motif_collection.getOverlap(subpeak_locus,'both'):
                out_degree_delta_dict[tf_name].append(signal_dict[subpeak_locus.__str__()])
        
    out_degree_table = [['TF_NAME','DELTA_MEAN','DELTA_MEDIAN','DELTA_STD']]
    for tf_name in tf_list:
        delta_line = [tf_name] + [numpy.mean(out_degree_delta_dict[tf_name]),numpy.median(out_degree_delta_dict[tf_name]),numpy.std(out_degree_delta_dict[tf_name])]
        out_degree_table.append(delta_line)


    utils.unParseTable(out_degree_table,output_path,'\t')
    print(output_path)
    return(output_path)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~CALCULATING CHANGES IN BRD4 OUT DEGREE BY TF CLUSTERS~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def tf_group_brd4_delta(crc_folder,chip_dataFile,analysis_name,tf_list,tf_group_name,y_brd4_list,o_brd4_list):

    '''
    calculates changes in brd4 out degree at each predicted motif occurrence
    '''
    crc_folder = utils.formatFolder(crc_folder,False)
    subpeak_bed_path = '%s%s_all_subpeak.bed' % (crc_folder,analysis_name)
    #direct the output to the crc folder
    signal_path = '%s%s_all_subpeak_signal.txt' % (crc_folder,analysis_name)
    
    #set final output
    output_path = '%s%s_all_subpeak_signal_%s.txt' % (crc_folder,analysis_name,tf_group_name)

    if len(tf_list) == 0:
        #next create a dictionary of all tfs in the system
        tf_table_path = '/storage/cylin/home/cl6/projects/NIBR_YvsO_cyl/crc/keratinocyte_combined/keratinocyte_combined_GENE_TF_TABLE.txt' 

        tf_table = utils.parseTable(tf_table_path,'\t')
        tf_list = utils.uniquify([line[0] for line in tf_table[1:]])
        tf_list.sort()

    #now bring in the signal table as a dictionary using the locus line as the id
    print('making log2 y vs o signal dict at subpeaks')
    signal_table = utils.parseTable(signal_path,'\t')
    signal_dict = defaultdict(float)
    
    #figure out columns for young and old
    o_columns = [signal_table[0].index(name) for name in o_brd4_list]
    y_columns = [signal_table[0].index(name) for name in y_brd4_list]
    for line in signal_table[1:]:
        o_signal = numpy.mean([float(line[col]) for col in o_columns])
        y_signal = numpy.mean([float(line[col]) for col in y_columns])
        delta = numpy.log2(y_signal/o_signal)
        signal_dict[line[1]] = delta


    #loading motifs
    print('making motif dictionaries')
    motif_dict = {}
    for tf_name in tf_list:
        print(tf_name)
        motif_bed_path = '%smotif_beds/%s_motifs.bed' % (crc_folder,tf_name)
        motif_bed = utils.parseTable(motif_bed_path,'\t')
        motif_loci = []
        for line in motif_bed[1:]:
            motif_locus = utils.Locus(line[0],line[1],line[2],'.')
            motif_loci.append(motif_locus)

        motif_collection = utils.LocusCollection(motif_loci)
        motif_dict[tf_name] = motif_collection


    #loading subpeak collection
    subpeak_table = utils.parseTable(subpeak_bed_path,'\t')
    #set the output table
    out_signal_table = [signal_table[0] + ['LENGTH'] + tf_list]
    for i in range(len(subpeak_table)):
        line = subpeak_table[i]
        line_locus = utils.Locus(line[0],int(line[1]),int(line[2]),'.')
        
        motif_line = [line_locus.len()]
        for tf_name in tf_list:
            motif_line.append(len(motif_dict[tf_name].getOverlap(line_locus)))
        out_signal_table.append([line[0]] + signal_table[(i+1)][1:] + motif_line)

    utils.unParseTable(out_signal_table,output_path,'\t')
    print(output_path)
    return(output_path)


        
    
#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
