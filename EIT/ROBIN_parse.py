# -*- coding: utf-8 -*-
import numpy as np
import base64
import os
import sys

# input_file : Holds name of file to parse
# mode = 0 -> Creates CSV's, Recommended for big files
# mode = 1 -> Fills in content in dictionary and ruturns it. Not recommended for big files
def ROBIN_parse(input_file, mode):

    # THIS IS WHAT NEEDS TO BE EDITED IF PAYLOADS CHANGES #    
    # SIGNAL DECLARATION ##########################################################
    ibis_bioz_blocksize = 256
    ibis_ecg_blocksize = 128
    ibis_accel_blocksize = 8
    ibis_freq_blocksize = 16
    ibis_ppg_blocksize = 25
    
    ibis_data_payload = np.dtype([('x','%dfloat32' %ibis_accel_blocksize),('y','%dfloat32' %ibis_accel_blocksize),('z','%dfloat32' %ibis_accel_blocksize),('ecg', '%dint32' %ibis_ecg_blocksize), ('eti','%dint32' %ibis_ecg_blocksize),('ecg2', '%dint32' %ibis_ecg_blocksize), ('Vi', '%dint32' %ibis_bioz_blocksize),('Vq', '%dint32' %ibis_bioz_blocksize),('freq', '%dint32' %ibis_freq_blocksize),('IR','%dint32' %ibis_ppg_blocksize),('red','%dint32' %ibis_ppg_blocksize)])
    
    ###############################################################################
    
    # Open file and read lines
    f = open(input_file, "rb")
    file_content = f.readlines()
    f.close()
    
    
    if mode == 0:    
        file_path = input_file[:-4]
        try:
            os.makedirs(file_path)
        except OSError:
            pass
        
        # Create csv files
        # Open up the CSV writer
        file_bioz = open(file_path + '/Bioz.csv', 'wb')
        file_freq = open(file_path + '/Freq.csv', 'wb')
        file_ecg = open(file_path + '/Ecg.csv', 'wb')
        file_accel = open(file_path + '/Accel.csv', 'wb')
        file_ppg = open(file_path + '/Ppg.csv', 'wb')
    
    # PARSE through the file
    num_payloads = len(file_content)-1
    
    # Create arrays to store data
    if mode == 1:    
        Data = dict()
        Data['vi'] = np.zeros(ibis_bioz_blocksize*num_payloads)
        Data['vq'] = np.zeros(ibis_bioz_blocksize*num_payloads)
        Data['ecg'] = np.zeros(ibis_ecg_blocksize*num_payloads)
        Data['eti'] = np.zeros(ibis_ecg_blocksize*num_payloads)
        Data['ecg2'] = np.zeros(ibis_ecg_blocksize*num_payloads)
        Data['x'] = np.zeros(ibis_accel_blocksize*num_payloads)
        Data['y'] = np.zeros(ibis_accel_blocksize*num_payloads)
        Data['z'] = np.zeros(ibis_accel_blocksize*num_payloads)
        Data['freq'] = np.zeros(ibis_freq_blocksize*num_payloads)
        Data['ir_ppg'] = np.zeros(ibis_ppg_blocksize*num_payloads)
        Data['red_ppg'] = np.zeros(ibis_ppg_blocksize*num_payloads)
    

    
    i = 0
    for line in file_content[1:]:   #Skip the first line because it's always old data
        #print 'Current payload read'
        try:
            current_payload = np.frombuffer(base64.b64decode(line[:-2]),dtype =ibis_data_payload)
            
            if mode == 1:
                Data['vi'][i*ibis_bioz_blocksize:(i+1)*ibis_bioz_blocksize] = current_payload['Vi']
                Data['vq'][i*ibis_bioz_blocksize:(i+1)*ibis_bioz_blocksize] = current_payload['Vq']
        
                Data['ecg'][i*ibis_ecg_blocksize:(i+1)*ibis_ecg_blocksize] = current_payload['ecg']
                Data['eti'][i*ibis_ecg_blocksize:(i+1)*ibis_ecg_blocksize] = current_payload['eti']
                Data['ecg2'][i*ibis_ecg_blocksize:(i+1)*ibis_ecg_blocksize] = current_payload['ecg2']
                
                Data['x'][i*ibis_accel_blocksize:(i+1)*ibis_accel_blocksize] = current_payload['x']
                Data['y'][i*ibis_accel_blocksize:(i+1)*ibis_accel_blocksize] = current_payload['y']
                Data['z'][i*ibis_accel_blocksize:(i+1)*ibis_accel_blocksize] = current_payload['z']
        
                Data['freq'][i*ibis_freq_blocksize:(i+1)*ibis_freq_blocksize] = current_payload['freq']
        
                Data['ir_ppg'][i*ibis_freq_blocksize:(i+1)*ibis_freq_blocksize] = current_payload['IR']
                Data['red_ppg'][i*ibis_freq_blocksize:(i+1)*ibis_freq_blocksize] = current_payload['red']
    
            if mode == 0:            
            # Write to CSV files
                BioZ = np.column_stack((current_payload['Vi'][0], current_payload['Vq'][0]))
                np.savetxt(file_bioz,BioZ)
                
                Ecg = np.column_stack((current_payload['ecg'][0], current_payload['eti'][0], current_payload['ecg2'][0]))
                np.savetxt(file_ecg,Ecg)
                #
                Accel = np.column_stack((current_payload['x'][0], current_payload['y'][0], current_payload['z'][0]))
                np.savetxt(file_accel,Accel)        
                
                np.savetxt(file_freq, current_payload['freq'][0])
                
                Ppg = np.column_stack((current_payload['IR'][0], current_payload['red'][0]))
                np.savetxt(file_ppg,Ppg)   
    
        except:
            pass 
        #Update progress bar    
        i = i + 1
    
    if mode == 0:
        file_bioz.close()
        file_ecg.close()
        file_accel.close()
        file_freq.close()
        file_ppg.close()
        
    
    if mode == 1:
        return Data


#Data = ROBIN_parse('ROBIN_10_ECG_and_BioZ_on_Seulki.TXT', 1)
#ROBIN_parse('ROBIN_10_ECG_and_BioZ_on_Seulki.TXT', 0)
argument = sys.argv
file_name = sys.argv[1]
ROBIN_parse(file_name, 0)