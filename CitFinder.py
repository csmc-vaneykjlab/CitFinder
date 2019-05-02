#!/usr/bin/env python

"""
==================================================================================
CitFinder: Identify & Validate Citrullinated Peptides using Mass Spectrometry Data
==================================================================================
Copyright (c) 2019 Ruining Liu & Vidya Venkatraman
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

from __future__ import print_function
import sys
import os
import csv
import getopt
import re
from numpy import *
import msproteomicstoolslib.format.speclib_db_lib as speclib_db_lib
  
#the following part is for getting modification site and 10 amino acid
def find(seq,database):#put \ before [ and ] this is for find the peptide (which remove the modification) in the database
    seq = seq.replace("[","\[")
    seq = seq.replace("]","\]")
    start_index_list = []
    for m in re.finditer(seq, database):
        start_index_list.append(m.start())
    return start_index_list 

def findMod(seq): #this is for find all the modification position in the peptide
    result_dict = {}
    i = 0
    index = 0
    for letter in seq:
        if letter == "[":
            start_point_site = i
            start_point = index
        if letter == "]":
            end_point = index
            amino_acid = seq[start_point-1:end_point+1]
            if result_dict.has_key(amino_acid):
                result_dict[amino_acid].append(start_point_site-1)
            else:
                result_dict[amino_acid] = [start_point_site-1]
            i = i -(end_point - start_point) - 1
        i += 1
        index += 1
    return result_dict

def getFasta(fasta):
    with open(fasta) as database:
        database_dict = {}
        proteinID = ''
        seq = ''
        for line in database:
            if line.startswith(">"):
                if "DECOY" in line:
                    if proteinID:
                        if database_dict.has_key(proteinID):
                            database_dict[proteinID].append(seq)
                        else:
                            database_dict[proteinID] = [seq]
                    seq = ''
                    proteinID = ''    
                    continue
                    
                else:
                    if proteinID:
                        if database_dict.has_key(proteinID):
                            database_dict[proteinID].append(seq)
                        else:
                            database_dict[proteinID] = [seq]
                        proteinID = line.strip(">\n").split(" ")[0]
                    else:
                        proteinID = line.strip(">\n").split(" ")[0]
                    seq = ''
            else:
                seq += line.strip("\n")
        if proteinID:
            if database_dict.has_key(proteinID):
                database_dict[proteinID].append(seq)
            else:
                database_dict[proteinID] = [seq]
            proteinID = line.strip(">\n").split(" ")[0]
        return database_dict

def getModSite(database_dict,proteinID,peptide,include_list):
    #include_list = "R[157]"#this is for the modification site and could be changed
    modInSeq_dict = {}
    modInSeq_all_dict = {}
    modification_site_list = []
    tenAA_list = []
    modInSeq_all_dict = findMod(peptide)
    for key,value in modInSeq_all_dict.iteritems():         
        if include_list in key:
            if modInSeq_dict.has_key(include_list):
                modInSeq_dict[include_list] += value
            else:
                modInSeq_dict[include_list] = value              
    peptide = re.sub('\[.*?\]', '', peptide)  #remove modification from peptide sequence   
    seqInDatabase_list  = find(peptide,database_dict[proteinID][0])
    if not seqInDatabase_list:
        return "",""
    if modInSeq_dict.has_key(include_list):
        for v_site in modInSeq_dict[include_list]:
            modi_site_amino_acid_list = []
            for seq_site in seqInDatabase_list:
                modification_site = int(seq_site) + int(v_site) +1
                start = modification_site-11
                end = modification_site+10
                start_flag = False
                end_flag = False
                seq_value = database_dict[proteinID][0]
                if start < 0:
                    start_flag = True
                if end > len(seq_value):
                    end_flag = True
                    end_difference = end - len(seq_value)
                if start_flag and not end_flag:
                    tenAA = abs(start)*'*' + seq_value[:end] 
                elif not start_flag and end_flag:
                    tenAA = seq_value[start:] + end_difference*'*'
                elif start_flag and end_flag:
                    tenAA = abs(start)*'*' + seq_value + end_difference*'*'
                else: 
                    tenAA = seq_value[start:end]
                modi_site_amino_acid_list.append(str(modification_site))
                tenAA_list.append(str(tenAA))
            modi_site_amino_acid = modi_site_amino_acid_list[0]
            modification_site_list.append(include_list[:include_list.find('[')]+modi_site_amino_acid)
    return ",".join(modification_site_list),",".join(tenAA_list)

def neutralLoss(sptxtfile):
    library_key = 99
    spectrastlib = speclib_db_lib.Library(library_key)
    offset = spectrastlib.get_first_offset(sptxtfile)
    last_offset = -100    
    sequence_dict = {}
    while ( offset - last_offset > 10) :
        last_offset = offset
        offset , spectrum = spectrastlib.read_sptxt_with_offset(sptxtfile,offset)
        sequence = spectrum.name.split('/')[0]
        charge = spectrum.name.split('/')[1]
        if 'R[157]' in sequence:
            site_specific_ion_b = re.sub('\[.*?\]', '', sequence).index('R')
            site_specific_ion_y = re.sub('\[.*?\]', '', sequence)[::-1].index('R')
            ion_list = []
            ion_43_list = []
            annotation_dict = {}#key:annotation,value:[peak,intensity]
            peaks = spectrum.get_peaks()        
            for peak in peaks :
                annotation = peak.peak_annotation
                for anno in annotation.split(','):
                    ion_part = anno.split('/')[0]
                    ion_list.append(ion_part)               
                    if '-43' in ion_part and ion_part[0] in ['b','y']:
                        site_ion = int(ion_part[1:ion_part.index('-')])
                        if (ion_part[0] == 'b' and site_ion > site_specific_ion_b) or (ion_part[0] == 'y' and site_ion > site_specific_ion_y):
                            ion_43_list.append(ion_part)
                    if annotation_dict.has_key(ion_part):
                        if float(peak.intensity) > float(annotation_dict[ion_part]['intensity']):
                            annotation_dict[ion_part]['peak'] = str(peak.peak)
                            annotation_dict[ion_part]['intensity'] = str(peak.intensity)
                    else:
                        annotation_dict[ion_part] = {'peak':str(peak.peak),'intensity':str(peak.intensity)}
            for ion_43 in set(ion_43_list):
                ion = ion_43.replace('-43','')
                if ion in ion_list: 
                    if sequence_dict.has_key(sequence):
                        ion_abs = ion_43.split(':')[0].split('i')[0].split('^')[0]
                        if sequence_dict[sequence].has_key(ion_abs):
                            if float(annotation_dict[ion_43]['intensity']) > float(sequence_dict[sequence][ion_abs][2]):
                                sequence_dict[sequence][ion_abs] = [ion,ion_43,annotation_dict[ion]['intensity'],annotation_dict[ion_43]['intensity']]   
                        else:
                            sequence_dict[sequence][ion_abs] = [ion,ion_43,annotation_dict[ion]['intensity'],annotation_dict[ion_43]['intensity']]
                             
                    else:
                        sequence_dict[sequence] = {}
                        ion_abs = ion_43.split(':')[0].split('i')[0].split('^')[0]
                        sequence_dict[sequence][ion_abs] = [ion,ion_43,annotation_dict[ion]['intensity'],annotation_dict[ion_43]['intensity']]
    return sequence_dict        
    
def lmedian(valarr):
  vals = sorted(valarr)
  if len(vals) % 2 == 1:
    return vals[(len(vals) + 1) // 2 - 1]
  else:
    return vals[len(vals) // 2 - 1]

def all_indices(value, qlist):
  indices = []
  idx = -1
  while True:
    try:
      idx = qlist.index(value, idx+1)
      indices.append(idx)
    except ValueError:
      break
  return indices
  
def readInput(file):
    try:
      sptxt_infile = open(file, 'r')
    except IOError:
      print(file, "not readable")
    peptide_dict_charge = {}#key:peptide key:rawspectrum key:rt  value: charge
    peptide_dict_charge_all = {}#key:peptide value: charge
    protein_dict = {}#key: peptide value: protein list
    rt_all_dict = {} #key: spectrum key: peptide value: rt list
    prob_dict = {} #key: spectrum key: peptide value: prob list
    peptide_charge_rt_dict = {} #key:peptide key:rawspectrum key:charge value:rt
    sptxt_header = []
    sptxt_block = []
    all_spectrum = []
    for sptxt_line in sptxt_infile:
      if sptxt_line[0] == "#":
        sptxt_header.append(sptxt_line)
      else:
        sptxt_block.append(sptxt_line)
        if sptxt_line == "\n":
          #peptide = sptxt_block[0].split("Name: ")[1].split("/")[0]
          peptide = sptxt_block[0].split("Name: ")[1].split("/")[0]
          libid = sptxt_block[1].split("LibID: ")[1].split("\n")[0]
          charge = sptxt_block[0].split("Name: ")[1].split("/")[1].split("\n")[0]
          mods = sptxt_block[6].split("Mods=")[1].split(" ")[0]
          protein = sptxt_block[6].split("Protein=")[1].split(" ")[0]
          spectrum = sptxt_block[6].split("RawSpectrum=")[1].split(".")[0]
          rt = float(sptxt_block[6].split("RetentionTime=")[1].split(",")[0])
          prob = float(sptxt_block[6].split("Prob=")[1].split(" ")[0])
          sptxt_block = []
          all_spectrum.append(spectrum)
          #gether peptide_dict_charge and peptide_dict_charge_all         
          if peptide_dict_charge_all.has_key(peptide):
              peptide_dict_charge_all[peptide].append(charge)
          else:
              peptide_dict_charge_all[peptide] = [charge]    
          if peptide_dict_charge.has_key(peptide):
              if peptide_dict_charge[peptide].has_key(spectrum):
                  if peptide_dict_charge[peptide][spectrum].has_key(rt):
                      peptide_dict_charge[peptide][spectrum][rt].append(charge)
                  else:
                      peptide_dict_charge[peptide][spectrum][rt] = [charge]                        
              else:
                  peptide_dict_charge[peptide][spectrum] = {}
                  peptide_dict_charge[peptide][spectrum][rt] = [charge]
          else:
              peptide_dict_charge[peptide] = {}
              peptide_dict_charge[peptide][spectrum] = {}
              peptide_dict_charge[peptide][spectrum][rt] = [charge]
          #gether protein_dict
          if protein_dict.has_key(peptide):
              protein_dict[peptide].append(protein) 
          else:
              protein_dict[peptide] = [protein] 
          #gether rt_all_dict
          if rt_all_dict.has_key(spectrum):
              if rt_all_dict[spectrum].has_key(peptide):
                  rt_all_dict[spectrum][peptide].append(rt)
              else:
                  rt_all_dict[spectrum][peptide] = [rt]     
          else:
              rt_all_dict[spectrum] = {}
              rt_all_dict[spectrum][peptide] = [rt]  
          #gether prob_dict
          if prob_dict.has_key(spectrum):
              if prob_dict[spectrum].has_key(peptide):
                  prob_dict[spectrum][peptide].append(prob)
              else:
                  prob_dict[spectrum][peptide] = [prob]     
          else:
              prob_dict[spectrum] = {}
              prob_dict[spectrum][peptide] = [prob]  
          #gether peptide_charge_rt_dict
          if peptide_charge_rt_dict.has_key(peptide):
              if peptide_charge_rt_dict[peptide].has_key(spectrum):
                  peptide_charge_rt_dict[peptide][spectrum][charge] = rt
              else:
                  peptide_charge_rt_dict[peptide][spectrum] = {}
                  peptide_charge_rt_dict[peptide][spectrum][charge] = rt
              
          else:
             peptide_charge_rt_dict[peptide] = {}
             peptide_charge_rt_dict[peptide][spectrum] = {}
             peptide_charge_rt_dict[peptide][spectrum][charge] = rt
    
    all_spectrum = list(set(all_spectrum))          
    sptxt_infile.close()  
    return peptide_charge_rt_dict,peptide_dict_charge, peptide_dict_charge_all, protein_dict, rt_all_dict, prob_dict, all_spectrum

def transferMedian(rt_all_dict,prob_dict):  
    rt_dict = {} #key: rawspectrum key: peptide value: Imedian rt
    rt_run_dict = {} #key: rawspectrum key: peptide value: sorted rt list?
    for rawspectrum in rt_all_dict:
        rt_dict[rawspectrum] = {}
        rt_run_dict[rawspectrum] = {}
        for peptide in rt_all_dict[rawspectrum]:
            rt = []
            for idx in all_indices(sorted(prob_dict[rawspectrum][peptide], reverse=True)[0],prob_dict[rawspectrum][peptide]):
                rt.append(rt_all_dict[rawspectrum][peptide][idx])
            rt_dict[rawspectrum][peptide] = lmedian(rt)
            rt_run_dict[rawspectrum][peptide] = rt   
    return rt_dict, rt_run_dict 
    
def readSkyline(file):
    try:
      skyline_file = open(file, 'r')
    except IOError:
      print(file, "not readable")
    replace_dict = {"R[+1]":"R[157]",'M[+16]':'M[147]','N[+1]':'N[115]','Q[+10]':'Q[129]','C[+57]':'C[160]'}
    skyline_report = csv.reader(skyline_file, delimiter=',')
    header_skyline = skyline_report.next()
    skyline_report_dict = {} #key: peptide, key:spectrum key: charge value: whole row
    mod_peptide_skyline_index = header_skyline.index('Peptide Modified Sequence')
    file_skyline_index = header_skyline.index('Replicate Name')
    charge_skyline_index = header_skyline.index('Precursor Charge')
    for row in skyline_report:
        peptide = row[mod_peptide_skyline_index]
        for old,new in replace_dict.iteritems():
            peptide = peptide.replace(old,new)
        file_name = row[file_skyline_index]
        charge = row[charge_skyline_index]
        if skyline_report_dict.has_key(peptide):
            if skyline_report_dict[peptide].has_key(file_name):
                skyline_report_dict[peptide][file_name][charge] = row
            else:
                skyline_report_dict[peptide][file_name] = {}
                skyline_report_dict[peptide][file_name][charge] = row              
        else:
            skyline_report_dict[peptide] = {}
            skyline_report_dict[peptide][file_name] = {}
            skyline_report_dict[peptide][file_name][charge] = row
    return skyline_report_dict,header_skyline

def skylineGoodCheck(skyline_row,header_skyline,mod_rt):
    flag_good = False
    idotp_index = header_skyline.index('Isotope Dot Product')
    rt_index = header_skyline.index('Best Retention Time')
    total_area_index = header_skyline.index('Total Area')
    fwhm_index = header_skyline.index('Max Fwhm')        
    average_mass_error_PPM_index = header_skyline.index('Average Mass Error PPM')
    if float(skyline_row[idotp_index]) >= 0.9 and abs(float(mod_rt)-float(skyline_row[rt_index])) <= 0.2 and float(skyline_row[total_area_index]) >= 100000 and float(skyline_row[total_area_index])/float(skyline_row[fwhm_index]) >=1000000 and abs(float(skyline_row[average_mass_error_PPM_index])) <= 5:
        flag_good = True
    return flag_good
    
def skylineOkayCheck(skyline_row,header_skyline,mod_rt):
    flag_okay = False
    score = 0
    idotp_index = header_skyline.index('Isotope Dot Product')
    rt_index = header_skyline.index('Best Retention Time')
    total_area_index = header_skyline.index('Total Area')
    fwhm_index = header_skyline.index('Max Fwhm')        
    average_mass_error_PPM_index = header_skyline.index('Average Mass Error PPM')
    if skyline_row[idotp_index] != '#N/A' and float(skyline_row[idotp_index]) >= 0.7:
        score += 1
    if skyline_row[idotp_index] != '#N/A' and float(skyline_row[idotp_index]) < 0.7:
        score = -100
    if skyline_row[rt_index] != '#N/A' and abs(float(mod_rt)-float(skyline_row[rt_index])) >= 1:
        score = -100
    if skyline_row[rt_index] != '#N/A' and abs(float(mod_rt)-float(skyline_row[rt_index])) <= 0.4: 
        score += 1
    if skyline_row[total_area_index] != '#N/A' and float(skyline_row[total_area_index]) >= 10000:
        score += 1
    if skyline_row[fwhm_index] != 0 and skyline_row[fwhm_index] != '#N/A' and float(skyline_row[total_area_index])/float(skyline_row[fwhm_index]) >=100000:
        score += 1
    if skyline_row[average_mass_error_PPM_index] != '#N/A' and abs(float(skyline_row[average_mass_error_PPM_index])) <= 10:
        score += 1
    if score >= 4:
        flag_okay = True
    return flag_okay

def skylineValidation(skyline_report_dict,peptide_charge_rt_dict,peptide_dict_charge,header_skyline,file_ori,mod_pep,charge_mod,mod_rt,unmod_rt,rt_shift):
    validation = ""   
    idotp_index = header_skyline.index('Isotope Dot Product')
    rt_index = header_skyline.index('Best Retention Time')
    total_area_index = header_skyline.index('Total Area')
    fwhm_index = header_skyline.index('Max Fwhm') 
    average_mass_error_PPM_index = header_skyline.index('Average Mass Error PPM')
    file_skyline_index = header_skyline.index('Replicate Name')
    file_ori_new = file_ori
    charge_mod_new = charge_mod
    mod_rt_new = mod_rt
    rt_shift_new = rt_shift
    rt_drift = 0
    total_area_fwhm = 0
    skyline_info = []
    skyline_row = {}
    if skyline_report_dict.has_key(mod_pep):
        if skyline_report_dict[mod_pep].has_key(file_ori):
            if skyline_report_dict[mod_pep][file_ori].has_key(charge_mod):
                skyline_row = skyline_report_dict[mod_pep][file_ori][charge_mod]
                if skyline_row[rt_index] != '#N/A' and skyline_row[idotp_index] != '#N/A' and skyline_row[total_area_index] != '#N/A' and skyline_row[fwhm_index] != 0 and skyline_row[fwhm_index] != '#N/A' and skyline_row[average_mass_error_PPM_index] != '#N/A':
                    if skylineGoodCheck(skyline_row,header_skyline,mod_rt):
                        validation = "Good" 
                    else:    
                        if skylineOkayCheck(skyline_row,header_skyline,mod_rt):
                            validation = 'Okay'
                        else:
                            validation = 'Bad'        
                else:                                  
                    if skylineOkayCheck(skyline_row,header_skyline,mod_rt):
                        validation = 'Okay'
                    else:
                        validation = 'Bad'
            else:
                validation = 'Bad'
        else:
            validation = 'Bad'
        
        if validation == 'Bad' and skyline_report_dict.has_key(mod_pep): #check other files
            all_values = skyline_report_dict[mod_pep].keys()#all spectrum
            file_rt = peptide_charge_rt_dict[mod_pep]
            all_spectrum = file_rt.keys()
            new_validation_yes_spectrum = ""
            idotp_yes = 0
            new_validation_soso_spectrum = ""
            idotp_soso = 0
            for spectrum_skyline in all_spectrum:
                if skyline_report_dict[mod_pep].has_key(spectrum_skyline):
                    value = skyline_report_dict[mod_pep][spectrum_skyline].values()[0]
                else:
                    continue
                charge_pep = file_rt[spectrum_skyline].keys()[0]
                new_mod_rt = round(float(file_rt[spectrum_skyline][charge_pep])/60,2)
                if new_mod_rt - unmod_rt >= 5:
                    if value[idotp_index] != '#N/A' and value[rt_index] != '#N/A' and value[total_area_index] != '#N/A' and value[fwhm_index] != '0' and value[fwhm_index] != '#N/A' and value[average_mass_error_PPM_index] != '#N/A':
                        if skylineGoodCheck(value,header_skyline,new_mod_rt):
                            if idotp_yes == 0:
                                idotp_yes = float(value[idotp_index]) 
                                new_validation_yes_spectrum = value
                            elif float(value[idotp_index]) > idotp_yes:
                                idotp_yes = float(value[idotp_index]) 
                                new_validation_yes_spectrum = value
                        else:
                            if skylineOkayCheck(value,header_skyline,new_mod_rt):
                                if idotp_soso == 0:
                                    idotp_soso = float(value[idotp_index]) 
                                    new_validation_soso_spectrum = value
                                elif float(value[idotp_index]) > idotp_soso:
                                    idotp_soso = float(value[idotp_index]) 
                                    new_validation_soso_spectrum = value
                    else:
                        if skylineOkayCheck(value,header_skyline,new_mod_rt):
                            if idotp_soso == 0:
                                idotp_soso = float(value[idotp_index]) 
                                new_validation_soso_spectrum = value
                            elif float(value[idotp_index]) > idotp_soso:
                                idotp_soso = float(value[idotp_index]) 
                                new_validation_soso_spectrum = value
            if new_validation_yes_spectrum:
                file_ori_new = '*'+ new_validation_yes_spectrum[file_skyline_index]                 
                mod_rt_new = file_rt[new_validation_yes_spectrum[file_skyline_index]].values()[0]
                charge_mod_new = peptide_dict_charge[mod_pep][new_validation_yes_spectrum[file_skyline_index]][mod_rt_new][0]
                mod_rt_new = round(float(mod_rt_new)/60,2)
                rt_shift_new = mod_rt_new - unmod_rt
                if new_validation_yes_spectrum[rt_index] != '#N/A':
                    rt_drift = float(mod_rt_new) - float(new_validation_yes_spectrum[rt_index]) 
                else:
                    rt_drift = '#N/A'
                if new_validation_yes_spectrum[total_area_index] != '#N/A' and new_validation_yes_spectrum[fwhm_index] != '#N/A':
                    total_area_fwhm =  float(new_validation_yes_spectrum[total_area_index])/float(new_validation_yes_spectrum[fwhm_index])
                else:
                    total_area_fwhm = '#N/A'
                skyline_info = [new_validation_yes_spectrum[idotp_index],new_validation_yes_spectrum[rt_index],rt_drift,new_validation_yes_spectrum[total_area_index],new_validation_yes_spectrum[fwhm_index],total_area_fwhm,new_validation_yes_spectrum[average_mass_error_PPM_index],'Good']
                return mod_rt_new,rt_shift_new,file_ori_new,charge_mod_new,skyline_info
            elif new_validation_soso_spectrum:
                file_ori_new = '*'+ new_validation_soso_spectrum[file_skyline_index]                 
                mod_rt_new = file_rt[new_validation_soso_spectrum[file_skyline_index]].values()[0]
                charge_mod_new = peptide_dict_charge[mod_pep][new_validation_soso_spectrum[file_skyline_index]][mod_rt_new][0]
                mod_rt_new = round(float(mod_rt_new)/60,2)
                rt_shift_new = mod_rt_new - unmod_rt
                if new_validation_soso_spectrum[rt_index] != '#N/A':
                    rt_drift = float(mod_rt_new) - float(new_validation_soso_spectrum[rt_index]) 
                else:
                    rt_drift = '#N/A'
                if new_validation_soso_spectrum[total_area_index] != '#N/A' and new_validation_soso_spectrum[fwhm_index] != '#N/A':
                    total_area_fwhm =  float(new_validation_soso_spectrum[total_area_index])/float(new_validation_soso_spectrum[fwhm_index])
                else:
                    total_area_fwhm = '#N/A'
                skyline_info = [new_validation_soso_spectrum[idotp_index],new_validation_soso_spectrum[rt_index],rt_drift,new_validation_soso_spectrum[total_area_index],new_validation_soso_spectrum[fwhm_index],total_area_fwhm,new_validation_soso_spectrum[average_mass_error_PPM_index],'Okay']
                return mod_rt_new,rt_shift_new,file_ori_new,charge_mod_new,skyline_info
            else:
                if skyline_row:
                    if skyline_row[rt_index] != '#N/A':
                        rt_drift = float(mod_rt_new) - float(skyline_row[rt_index]) 
                    else:
                        rt_drift = '#N/A'
                    if skyline_row[total_area_index] != '#N/A' and skyline_row[fwhm_index] != '#N/A':
                        total_area_fwhm =  float(skyline_row[total_area_index])/float(skyline_row[fwhm_index])
                    else:
                        total_area_fwhm = '#N/A'
                    skyline_info = [skyline_row[idotp_index],skyline_row[rt_index],rt_drift,skyline_row[total_area_index],skyline_row[fwhm_index],total_area_fwhm,skyline_row[average_mass_error_PPM_index],'Bad']
                    return mod_rt_new,rt_shift_new,file_ori_new,charge_mod_new,skyline_info    
                else:
                    skyline_info = ['#N/A','#N/A','#N/A','#N/A','#N/A','#N/A','#N/A','Bad']
                    return mod_rt_new,rt_shift_new,file_ori_new,charge_mod_new,skyline_info                    
        else:
            if skyline_row[rt_index] != '#N/A':
                rt_drift = float(mod_rt_new) - float(skyline_row[rt_index]) 
            else:
                rt_drift = '#N/A'
            if skyline_row[total_area_index] != '#N/A' and skyline_row[fwhm_index] != '#N/A':
                total_area_fwhm =  float(skyline_row[total_area_index])/float(skyline_row[fwhm_index])
            else:
                total_area_fwhm = '#N/A'
            skyline_info = [skyline_row[idotp_index],skyline_row[rt_index],rt_drift,skyline_row[total_area_index],skyline_row[fwhm_index],total_area_fwhm,skyline_row[average_mass_error_PPM_index],validation]
            return mod_rt_new,rt_shift_new,file_ori_new,charge_mod_new,skyline_info
    else:
        skyline_info = ['#N/A','#N/A','#N/A','#N/A','#N/A','#N/A','#N/A','Bad']
        return mod_rt_new,rt_shift_new,file_ori_new,charge_mod_new,skyline_info

def validationReport(rt_all_dict,peptide_dict,peptide_charge_rt_dict,peptide_dict_charge,peptide_dict_charge_all,rt_dict,rt_run_dict,protein_dict,database_dict,sample_type,file_pair,include_list,sequence_dict,skyline_report_dict,header_skyline,rtShift_flag,format_peptide,format_site):
  try:
    filename, ext = os.path.splitext(file_pair)
    pair_header = ['Protein','Mod Peptide','Modification Site','10 Amino Acid','Peptide Pair','Mod RT','Unmod RT','RT Shift','Mod Charge State','Unmod Charge State','Charge Shift','Total Neutral Loss Ions ','Neutral Loss Pairs','Other Mods RT','Mod-Unmod File Pairs'] + sample_type
    if skyline_report_dict:
        pair_header += ['Isotope Dot Product','Best Retention Time','Retention Time Drift','Total Area','Max Fwhm','Total Area/Fwhm','Average Mass Error PPM','Skyline Validation']
    if format_peptide:
        pair_peptide = csv.writer(open(filename+"_peptide"+ext,'w'),delimiter=',')   
        pair_peptide.writerow(pair_header) 
    if format_site:
        pair_site = csv.writer(open(filename+"_site"+ext,'w'),delimiter=',') 
        pair_site.writerow(pair_header)
    
    ind = {}
    peptide_list = []
    peptide_r_list = []
    peptide_dict_pep = {}
    for spectrum in rt_dict:
      for peptide in rt_dict[spectrum]:
        if peptide not in ind:
          ind[peptide] = []
        ind[peptide].append(spectrum)
    for peptide in ind:
      peptide_list.append(peptide)
      if include_list in peptide:
          peptide_r_list.append(peptide) 
      key = re.sub('\[.*?\]', '', peptide)   
      if peptide_dict_pep.has_key(key):
          peptide_dict_pep[key].append(peptide)
      else:
          peptide_dict_pep[key] = [peptide]
      rt = []
      rt_run_all = []
      rt_run_median = []
      rt_run_mean = []
      rt_run_sd = []
      spectrum_list = [] 
      protein = protein_dict[peptide]
      for spectrum in ind[peptide]:
        spectrum_list.append(spectrum)
        rt.append(rt_dict[spectrum][peptide])
        for rt_value in rt_all_dict[spectrum][peptide]:
            rt_run_all.append(spectrum+":"+str(rt_value))
        rt_run_median.append(spectrum+":"+str(lmedian(rt_run_dict[spectrum][peptide])))        
        rt_run_mean.append(spectrum+":"+str(mean(rt_run_dict[spectrum][peptide])))
        rt_run_sd.append(spectrum+":"+str(std(rt_run_dict[spectrum][peptide])))
      protein_list = protein[0].split("/")
      for i in range(1,len(protein_list)):
          modificationSite,tenaa = getModSite(database_dict,protein_list[i],peptide,include_list)
          if modificationSite:
              break
      peptide_dict[peptide] = {}
      peptide_dict[peptide]["protein"] = protein[0]
      peptide_dict[peptide]["rt"] = lmedian(rt)
      peptide_dict[peptide]["spectrum"] = ";".join(rt_run_all)
      peptide_dict[peptide]["modificationSite"] = modificationSite
      peptide_dict[peptide]["tenaa"] = tenaa
      for spectrum_info in rt_run_median:
          if str(lmedian(rt)) in spectrum_info:               
              peptide_dict[peptide]['median_spectrum'] = spectrum_info.split(':')[0] 
    for peptide_r in peptide_r_list:
        if sequence_dict.has_key(peptide_r):
            neutral_loss_res = []
            v = sequence_dict[peptide_r]
            max_intensity =  max(v.keys())
            neutral_loss_count = len(v)
            for key,value in v.iteritems():
                neutral_loss_res.append(':'.join([value[0]+'('+value[2]+')',value[1]+'('+value[3]+')']))
        else:
            neutral_loss_count = ""
            neutral_loss_res = []
        protein_value = peptide_dict[peptide_r]["protein"]
        Modification_site_value = peptide_dict[peptide_r]["modificationSite"]
        aa_value = peptide_dict[peptide_r]["tenaa"]
        Modification_site_value_site = Modification_site_value.split(',')
        aa_value_site = aa_value.split(',')
        charge_mod_unique = peptide_dict_charge[peptide_r][peptide_dict[peptide_r]['median_spectrum']][peptide_dict[peptide_r]['rt']][0]
        charge_mod = peptide_dict_charge_all[peptide_r]
        sample_num = []
        for sample in sample_type:
            sample_num.append(peptide_dict[peptide_r]["spectrum"].count(sample))                 
        unmod = peptide_r.replace(include_list,include_list[0])#it can only target to one modification. For example, only R[157] or only N[115]
        naked = re.sub('\[.*?\]', '', peptide_r)
        rt_mod = round(float(peptide_dict[peptide_r]["rt"])/60,2)
        if unmod in peptide_list: #it means it has a pair
            rt_unmod = round(float(peptide_dict[unmod]["rt"])/60,2)
            rt_shift = rt_mod - rt_unmod 
            if rtShift_flag:
                if rt_shift < 5:
                    continue
            charge_unmod_unique = peptide_dict_charge[unmod][peptide_dict[unmod]['median_spectrum']][peptide_dict[unmod]['rt']][0]
            charge_unmod = peptide_dict_charge_all[unmod]
            if set(charge_mod) == set(charge_unmod):
                charge_shift_value = "Match"
            else:
                result_charge_shift = []
                for chargeMod in charge_mod:
                    for chargeUnmod in charge_unmod:
                        if chargeMod < chargeUnmod:
                            result_charge_shift.append("Loss")
                        elif chargeMod == chargeUnmod:
                            result_charge_shift.append("Match")
                charge_shift_value = ",".join(set(result_charge_shift))
            if skyline_report_dict:
                if rt_shift >= 5:
                    mod_rt_new,rt_shift_new,file_ori_new,charge_mod_new,skyline_info = skylineValidation(skyline_report_dict,peptide_charge_rt_dict,peptide_dict_charge,header_skyline,peptide_dict[peptide_r]['median_spectrum'],peptide_r,charge_mod_unique,rt_mod,rt_unmod,rt_shift)
                    if format_peptide:
                        row_peptide = [protein_value,peptide_r,Modification_site_value,aa_value,unmod,str(mod_rt_new),str(rt_unmod),str(rt_shift_new),','.join(set(charge_mod)),','.join(set(charge_unmod)),charge_shift_value,neutral_loss_count,';'.join(neutral_loss_res),'',file_ori_new+': '+str(charge_mod_new)+', '+peptide_dict[unmod]['median_spectrum']+': '+str(charge_unmod_unique)] + sample_num + skyline_info
                    if format_site:
                        for i in range(len(Modification_site_value_site)):
                            row_site = [protein_value,peptide_r,Modification_site_value_site[i],aa_value_site[i],unmod,str(mod_rt_new),str(rt_unmod),str(rt_shift_new),','.join(set(charge_mod)),','.join(set(charge_unmod)),charge_shift_value,neutral_loss_count,';'.join(neutral_loss_res),'',file_ori_new+': '+str(charge_mod_new)+', '+peptide_dict[unmod]['median_spectrum']+': '+str(charge_unmod_unique)] + sample_num + skyline_info
                            pair_site.writerow(row_site)
                    
                else:
                    if format_peptide:
                        row_peptide = [protein_value,peptide_r,Modification_site_value,aa_value,unmod,str(rt_mod),str(rt_unmod),str(rt_shift),','.join(set(charge_mod)),','.join(set(charge_unmod)),charge_shift_value,neutral_loss_count,';'.join(neutral_loss_res),'',peptide_dict[peptide_r]['median_spectrum']+': '+str(charge_mod_unique)+', '+peptide_dict[unmod]['median_spectrum']+': '+str(charge_unmod_unique)] + sample_num + ['','','','','','','','#N/A']
                    if format_site:
                        for i in range(len(Modification_site_value_site)):
                            row_site = [protein_value,peptide_r,Modification_site_value_site[i],aa_value_site[i],unmod,str(rt_mod),str(rt_unmod),str(rt_shift),','.join(set(charge_mod)),','.join(set(charge_unmod)),charge_shift_value,neutral_loss_count,';'.join(neutral_loss_res),'',peptide_dict[peptide_r]['median_spectrum']+': '+str(charge_mod_unique)+', '+peptide_dict[unmod]['median_spectrum']+': '+str(charge_unmod_unique)] + sample_num + ['','','','','','','','#N/A']
                            pair_site.writerow(row_site)
            else:
                if format_peptide:
                    row_peptide = [protein_value,peptide_r,Modification_site_value,aa_value,unmod,str(rt_mod),str(rt_unmod),str(rt_shift),','.join(set(charge_mod)),','.join(set(charge_unmod)),charge_shift_value,neutral_loss_count,';'.join(neutral_loss_res),'',peptide_dict[peptide_r]['median_spectrum']+': '+str(charge_mod_unique)+', '+peptide_dict[unmod]['median_spectrum']+': '+str(charge_unmod_unique)] + sample_num
                if format_site:
                    for i in range(len(Modification_site_value_site)):
                        row_site = [protein_value,peptide_r,Modification_site_value_site[i],aa_value_site[i],unmod,str(rt_mod),str(rt_unmod),str(rt_shift),','.join(set(charge_mod)),','.join(set(charge_unmod)),charge_shift_value,neutral_loss_count,';'.join(neutral_loss_res),'',peptide_dict[peptide_r]['median_spectrum']+': '+str(charge_mod_unique)+', '+peptide_dict[unmod]['median_spectrum']+': '+str(charge_unmod_unique)] + sample_num
                        pair_site.writerow(row_site)
            if format_peptide:
                pair_peptide.writerow(row_peptide)             
        else:
            if len(findMod(peptide_r).keys()) == 1 : #this is the case we will bring it back
                count_mod = peptide_r.count(include_list)
                other_peptide_pair = {}
                for all_peptide in peptide_dict_pep[naked]:
                    if all_peptide.count(include_list) < count_mod and len(findMod(all_peptide).keys()) == 1:
                        rt_key = round(float(peptide_dict[all_peptide]["rt"])/60,2)
                        other_peptide_pair[rt_key] = all_peptide
                if other_peptide_pair:
                    all_peptide_selected_list = [max(other_peptide_pair.values(), key=len)]
                    len_selected = len(all_peptide_selected_list[0])
                    for peptide in other_peptide_pair.itervalues():
                        if len(peptide) == len_selected and peptide not in all_peptide_selected_list:
                            all_peptide_selected_list.append(peptide)  
                    if  len(all_peptide_selected_list) > 1:#multiple possible pairs
                        rt_list = []
                        for key,value in other_peptide_pair.iteritems():
                            if value in all_peptide_selected_list:
                                rt_list.append(key)
                        if max(rt_list) - min(rt_list) > 1: #the min and min has greater than 1 min rt shift, then kick it out
                            other_peptide = []
                            if rtShift_flag:                            
                                continue
                            for all_peptide in peptide_dict_pep[naked]:
                                if all_peptide == peptide_r:
                                    continue
                                else:
                                    rt_other = round(float(peptide_dict[all_peptide]["rt"])/60,2)
                                    other_peptide.append(all_peptide+'('+str(rt_other)+')')
                            if skyline_report_dict:
                                if format_peptide:
                                    row_peptide = [protein_value,peptide_r,Modification_site_value,aa_value,'',str(rt_mod),'','',','.join(set(charge_mod)),'','',neutral_loss_count,';'.join(neutral_loss_res),'; '.join(other_peptide),peptide_dict[peptide_r]['median_spectrum']] + sample_num + ['','','','','','','','#N/A']
                                if format_site:
                                    for i in range(len(Modification_site_value_site)):
                                        row_site = [protein_value,peptide_r,Modification_site_value_site[i],aa_value_site[i],'',str(rt_mod),'','',','.join(set(charge_mod)),'','',neutral_loss_count,';'.join(neutral_loss_res),'; '.join(other_peptide),peptide_dict[peptide_r]['median_spectrum']] + sample_num + ['','','','','','','','#N/A']
                                        pair_site.writerow(row_site)
                            else:
                                if format_peptide:
                                    row_peptide = [protein_value,peptide_r,Modification_site_value,aa_value,'',str(rt_mod),'','',','.join(set(charge_mod)),'','',neutral_loss_count,';'.join(neutral_loss_res),'; '.join(other_peptide),peptide_dict[peptide_r]['median_spectrum']] + sample_num
                                if format_site:
                                    for i in range(len(Modification_site_value_site)):
                                        row_site = [protein_value,peptide_r,Modification_site_value_site[i],aa_value_site[i],'',str(rt_mod),'','',','.join(set(charge_mod)),'','',neutral_loss_count,';'.join(neutral_loss_res),'; '.join(other_peptide),peptide_dict[peptide_r]['median_spectrum']] + sample_num
                                        pair_site.writerow(row_site)
                            if format_peptide:
                                pair_peptide.writerow(row_peptide)
                            continue
                        else:
                            rt_shift_list_check = []
                            for peptide in all_peptide_selected_list:
                                rt_peptide_unmod = round(float(peptide_dict[peptide]["rt"])/60,2)
                                if abs(rt_mod-rt_peptide_unmod) >= 5:
                                    rt_shift_list_check.append("shift") 
                                else:
                                    rt_shift_list_check.append("no shift")  
                            if len(set(rt_shift_list_check)) == 1:# will keep those ones
                                all_peptide_selected = all_peptide_selected_list[0]
                            else:#throw it away cause they are in different trend  
                                other_peptide = []
                                if rtShift_flag:
                                    continue
                                for all_peptide in peptide_dict_pep[naked]:
                                    if all_peptide == peptide_r:
                                        continue
                                    else:
                                        rt_other = round(float(peptide_dict[all_peptide]["rt"])/60,2)
                                        other_peptide.append(all_peptide+'('+str(rt_other)+')')
                                if skyline_report_dict:
                                    if format_peptide:
                                        row_peptide = [protein_value,peptide_r,Modification_site_value,aa_value,'',str(rt_mod),'','',','.join(set(charge_mod)),'','',neutral_loss_count,';'.join(neutral_loss_res),'; '.join(other_peptide),peptide_dict[peptide_r]['median_spectrum']] + sample_num + ['','','','','','','','#N/A']
                                    if format_site:
                                        for i in range(len(Modification_site_value_site)):
                                            row_site = [protein_value,peptide_r,Modification_site_value_site[i],aa_value_site[i],'',str(rt_mod),'','',','.join(set(charge_mod)),'','',neutral_loss_count,';'.join(neutral_loss_res),'; '.join(other_peptide),peptide_dict[peptide_r]['median_spectrum']] + sample_num + ['','','','','','','','#N/A']
                                            pair_site.writerow(row_site)
                                else:
                                    if format_peptide:
                                        row_peptide = [protein_value,peptide_r,Modification_site_value,aa_value,'',str(rt_mod),'','',','.join(set(charge_mod)),'','',neutral_loss_count,';'.join(neutral_loss_res),'; '.join(other_peptide),peptide_dict[peptide_r]['median_spectrum']] + sample_num 
                                    if format_site:
                                        for i in range(len(Modification_site_value_site)): 
                                            row_site = [protein_value,peptide_r,Modification_site_value_site[i],aa_value_site[i],'',str(rt_mod),'','',','.join(set(charge_mod)),'','',neutral_loss_count,';'.join(neutral_loss_res),'; '.join(other_peptide),peptide_dict[peptide_r]['median_spectrum']] + sample_num 
                                            pair_site.writerow(row_site)
                                if format_peptide:      
                                    pair_peptide.writerow(row_peptide)
                                continue  
                    else:
                        all_peptide_selected = all_peptide_selected_list[0]
                    rt_all_peptide = round(float(peptide_dict[all_peptide_selected]["rt"])/60,2)
                    charge_all_peptide_unique = peptide_dict_charge[all_peptide_selected][peptide_dict[all_peptide_selected]['median_spectrum']][peptide_dict[all_peptide_selected]['rt']][0]
                    charge_all_peptide = peptide_dict_charge_all[all_peptide_selected]
                    rt_shift = rt_mod - rt_all_peptide
                    if set(charge_mod) == set(charge_all_peptide):
                        charge_shift_value = "Match"
                    else:
                        result_charge_shift = []
                        for chargeMod in charge_mod:
                            for chargeUnmod in charge_all_peptide:
                                if chargeMod < chargeUnmod:
                                    result_charge_shift.append("Loss")
                                elif chargeMod == chargeUnmod:
                                    result_charge_shift.append("Match")
                        charge_shift_value = ",".join(set(result_charge_shift))
                    if rtShift_flag:
                        if rt_shift < 5:
                            continue
                    if skyline_report_dict:
                        if rt_shift >= 5:
                            mod_rt_new,rt_shift_new,file_ori_new,charge_mod_new,skyline_info = skylineValidation(skyline_report_dict,peptide_charge_rt_dict,peptide_dict_charge,header_skyline,peptide_dict[peptide_r]['median_spectrum'],peptide_r,charge_mod_unique,rt_mod,rt_all_peptide,rt_shift)
                            if format_peptide:
                                row_peptide = [protein_value,peptide_r,Modification_site_value,aa_value,all_peptide_selected,str(mod_rt_new),str(rt_all_peptide),str(rt_shift_new),','.join(set(charge_mod)),','.join(set(charge_all_peptide)),charge_shift_value,neutral_loss_count,';'.join(neutral_loss_res),'',file_ori_new+': '+str(charge_mod_new)+', '+peptide_dict[all_peptide_selected]['median_spectrum']+': '+str(charge_all_peptide_unique)] + sample_num + skyline_info
                            if format_site:
                                for i in range(len(Modification_site_value_site)):  
                                    row_site = [protein_value,peptide_r,Modification_site_value_site[i],aa_value_site[i],all_peptide_selected,str(mod_rt_new),str(rt_all_peptide),str(rt_shift_new),','.join(set(charge_mod)),','.join(set(charge_all_peptide)),charge_shift_value,neutral_loss_count,';'.join(neutral_loss_res),'',file_ori_new+': '+str(charge_mod_new)+', '+peptide_dict[all_peptide_selected]['median_spectrum']+': '+str(charge_all_peptide_unique)] + sample_num + skyline_info   
                                    pair_site.writerow(row_site)
                        else: 
                            if format_peptide:
                                row_peptide = [protein_value,peptide_r,Modification_site_value,aa_value,all_peptide_selected,str(rt_mod),str(rt_all_peptide),str(rt_shift),','.join(set(charge_mod)),','.join(set(charge_all_peptide)),charge_shift_value,neutral_loss_count,';'.join(neutral_loss_res),'',peptide_dict[peptide_r]['median_spectrum']+': '+str(charge_mod_unique)+', '+peptide_dict[all_peptide_selected]['median_spectrum']+': '+str(charge_all_peptide_unique)] + sample_num + ['','','','','','','','#N/A']
                            if format_site:
                                for i in range(len(Modification_site_value_site)):  
                                    row_site = [protein_value,peptide_r,Modification_site_value_site[i],aa_value_site[i],all_peptide_selected,str(rt_mod),str(rt_all_peptide),str(rt_shift),','.join(set(charge_mod)),','.join(set(charge_all_peptide)),charge_shift_value,neutral_loss_count,';'.join(neutral_loss_res),'',peptide_dict[peptide_r]['median_spectrum']+': '+str(charge_mod_unique)+', '+peptide_dict[all_peptide_selected]['median_spectrum']+': '+str(charge_all_peptide_unique)] + sample_num + ['','','','','','','','#N/A']
                                    pair_site.writerow(row_site)
                    else:
                        if format_peptide:
                            row_peptide = [protein_value,peptide_r,Modification_site_value,aa_value,all_peptide_selected,str(rt_mod),str(rt_all_peptide),str(rt_shift),','.join(set(charge_mod)),','.join(set(charge_all_peptide)),charge_shift_value,neutral_loss_count,';'.join(neutral_loss_res),'',peptide_dict[peptide_r]['median_spectrum']+': '+str(charge_mod_unique)+', '+peptide_dict[all_peptide_selected]['median_spectrum']+': '+str(charge_all_peptide_unique)] + sample_num       
                        if format_site:
                            for i in range(len(Modification_site_value_site)):         
                                row_site = [protein_value,peptide_r,Modification_site_value_site[i],aa_value_site[i],all_peptide_selected,str(rt_mod),str(rt_all_peptide),str(rt_shift),','.join(set(charge_mod)),','.join(set(charge_all_peptide)),charge_shift_value,neutral_loss_count,';'.join(neutral_loss_res),'',peptide_dict[peptide_r]['median_spectrum']+': '+str(charge_mod_unique)+', '+peptide_dict[all_peptide_selected]['median_spectrum']+': '+str(charge_all_peptide_unique)] + sample_num 
                                pair_site.writerow(row_site)
                    if format_peptide:  
                        pair_peptide.writerow(row_peptide)                                                  
                else:
                    other_peptide = []
                    if rtShift_flag:
                        continue
                    for all_peptide in peptide_dict_pep[naked]:
                        if all_peptide == peptide_r:
                            continue
                        else:
                            rt_other = round(float(peptide_dict[all_peptide]["rt"])/60,2)
                            other_peptide.append(all_peptide+'('+str(rt_other)+')')                            
                    if skyline_report_dict:
                        if format_peptide:
                            row_peptide = [protein_value,peptide_r,Modification_site_value,aa_value,'',str(rt_mod),'','',','.join(set(charge_mod)),'','',neutral_loss_count,';'.join(neutral_loss_res),'; '.join(other_peptide),peptide_dict[peptide_r]['median_spectrum']] + sample_num +  ['','','','','','','','#N/A']
                        if format_site:
                            for i in range(len(Modification_site_value_site)):  
                                row_site = [protein_value,peptide_r,Modification_site_value_site[i],aa_value_site[i],'',str(rt_mod),'','',','.join(set(charge_mod)),'','',neutral_loss_count,';'.join(neutral_loss_res),'; '.join(other_peptide),peptide_dict[peptide_r]['median_spectrum']] + sample_num +  ['','','','','','','','#N/A']
                                pair_site.writerow(row_site)
                    else: 
                        if format_peptide:   
                            row_peptide = [protein_value,peptide_r,Modification_site_value,aa_value,'',str(rt_mod),'','',','.join(set(charge_mod)),'','',neutral_loss_count,';'.join(neutral_loss_res),'; '.join(other_peptide),peptide_dict[peptide_r]['median_spectrum']] + sample_num    
                        if format_site:
                            for i in range(len(Modification_site_value_site)):
                                row_site = [protein_value,peptide_r,Modification_site_value_site[i],aa_value_site[i],'',str(rt_mod),'','',','.join(set(charge_mod)),'','',neutral_loss_count,';'.join(neutral_loss_res),'; '.join(other_peptide),peptide_dict[peptide_r]['median_spectrum']] + sample_num      
                                pair_site.writerow(row_site)   
                    if format_peptide:  
                        pair_peptide.writerow(row_peptide)
                
            else:
                other_peptide = []
                if rtShift_flag:
                    continue
                for all_peptide in peptide_dict_pep[naked]:
                    if all_peptide == peptide_r:
                        continue
                    else:
                        rt_other = round(float(peptide_dict[all_peptide]["rt"])/60,2)
                        other_peptide.append(all_peptide+'('+str(rt_other)+')')
                if skyline_report_dict:
                    if format_peptide:  
                        row_peptide = [protein_value,peptide_r,Modification_site_value,aa_value,'',str(rt_mod),'','',','.join(set(charge_mod)),'','',neutral_loss_count,';'.join(neutral_loss_res),'; '.join(other_peptide),peptide_dict[peptide_r]['median_spectrum']] + sample_num + ['','','','','','','','#N/A']
                    if format_site:
                        for i in range(len(Modification_site_value_site)):
                            row_site = [protein_value,peptide_r,Modification_site_value_site[i],aa_value_site[i],'',str(rt_mod),'','',','.join(set(charge_mod)),'','',neutral_loss_count,';'.join(neutral_loss_res),'; '.join(other_peptide),peptide_dict[peptide_r]['median_spectrum']] + sample_num + ['','','','','','','','#N/A']
                            pair_site.writerow(row_site)   
                else:   
                    if format_peptide:    
                        row_peptide = [protein_value,peptide_r,Modification_site_value,aa_value,'',str(rt_mod),'','',','.join(set(charge_mod)),'','',neutral_loss_count,';'.join(neutral_loss_res),'; '.join(other_peptide),peptide_dict[peptide_r]['median_spectrum']] + sample_num  
                    if format_site:
                        for i in range(len(Modification_site_value_site)):   
                            row_site = [protein_value,peptide_r,Modification_site_value_site[i],aa_value_site[i],'',str(rt_mod),'','',','.join(set(charge_mod)),'','',neutral_loss_count,';'.join(neutral_loss_res),'; '.join(other_peptide),peptide_dict[peptide_r]['median_spectrum']] + sample_num  
                if format_peptide:        
                    pair_peptide.writerow(row_peptide)                      
  except IOError:
    print(file, "not readable")          
            
# MAIN
def main(argv):
  splib_input = False
  rtShift_flag = True
  help = False
  fasta = ""
  all_spectrum = []
  sample_type = []
  include_list = []
  peptide_dict = {} #it has rt, protein, spectrum
  database_dict = {}
  skyline_report_dict = {}
  header_skyline = ""
  file_pair = ""
  format = []
  format_peptide = True
  format_site = False
  
  try:
    opts, args = getopt.getopt(argv, "i:o:f:g:m:s:r:F:",["in=","out=","fasta","grouping","modification","skyline","rtShift","format"])
  except getopt.GetoptError:
    help = True
    opts =  ( ("",""),)

  for opt, arg in opts:
    if opt in ("-i","--in"):
        splib_input = arg
    elif opt in ("-o","--out"):
        file_pair = arg
    elif opt in ("-f","--fasta"):
        fasta = arg
    elif opt in ("-g","--grouping"):
        sample_type= arg.split(",")
    elif opt in ("-m","--modification"):
        include_list = arg
        if '[' not in include_list:
            include_list = include_list + '['
    elif opt in ("-s","--skyline"):
        skyline_report = arg
        skyline_report_dict,header_skyline = readSkyline(skyline_report)
    elif opt in ("-r","--rtShift"):
        if arg in ['True','False']:
            if arg == 'False':
                rtShift_flag = False
        else:
            help = True 
    elif opt in ("-F","--format"):
        format = arg.split(",")
        format_peptide = False
        format_site = False
        for f in format:
            if f == "Peptide":
                format_peptide = True
            elif f == "Site":
                format_site = True
            else:
                help = True         
            
  if help or not splib_input or not file_pair or not fasta or not include_list:
    print("Cit Finder")
    print("---------------------------------------------------------------------------------------------")
    print("Usage:     CitFinder.py -i non_consensus_library.[splib/sptxt] -o validation_output.csv")
    print("Input:     SpectraST non_consensus_library.splib in txt format")
    print("Output:    Modified peptides pairs with RT information, neutral loss and skyline validation results.")
    print("Argument:  -i [--in]: input splib file")
    print("           -o [--out]: output csv file")
    print("           -f [--fasta]: specify the fasta file for modification site and 10 amino acid information")
    print("           -m [--modification]: specify the modification. Please specify one targed mod at a time. eg. R[157] OR R ONLY when both of the fasta file and modification are specified, pair file will be reported.")
    print("           (optional) -g [--grouping]: specify the grouping and comma seprate them. For example: Heart,Lung,Liver,Muscle,Kidney,Brain")
    print("           (optional) -s [--skyline]: Skyline report for validation")
    print("           (optional) -r [--rtShift]: If rt shift is True, it will only provide the modifed peptide pairs with >= 5 mins rt shift.Otherwise, it will provide all the modified peptide pairs. Default: True")
    print("           (optional) -F [--format]: The output format options are Peptide, Site. Please use comma to separate the multiple options. Default: Peptide")
    print("Important: The splib need to be in txt format!")
    print("Contact:   Ruining Liu <ruining.liu@cshs.org>")
    sys.exit()
  #get fasta database
  database_dict = getFasta(fasta)
  #read splib input
  peptide_charge_rt_dict,peptide_dict_charge, peptide_dict_charge_all, protein_dict, rt_all_dict, prob_dict, all_spectrum = readInput(splib_input)
  #get Median rt and probability
  rt_dict, rt_run_dict  = transferMedian(rt_all_dict,prob_dict) 
  #get neutral loss
  sequence_dict = neutralLoss(splib_input)
  #read skyline report
  if skyline_report_dict:
      skyline_report_dict,header_skyline = readSkyline(skyline_report)
  # write report
  if not sample_type:
      sample_type = sorted(all_spectrum)   
  validationReport(rt_all_dict,peptide_dict,peptide_charge_rt_dict,peptide_dict_charge,peptide_dict_charge_all,rt_dict,rt_run_dict,protein_dict,database_dict,sample_type,file_pair,include_list,sequence_dict,skyline_report_dict,header_skyline,rtShift_flag,format_peptide,format_site)
 
if __name__ == "__main__":
  main(sys.argv[1:])

