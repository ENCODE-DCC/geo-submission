import requests
import json
from urllib.parse import urlsplit, urlunsplit, parse_qs
import sys
from time import sleep, strftime
import pandas as pd
import time
import ast
import argparse
import os
from os.path import splitext


def main():
    t0 = time.process_time()
    args = get_args()
    # phase 1 - collect all experiments submitted so far.
    submitted_experiments = set()

    with open(args.infile) as exp_f:
        for l in exp_f:
            # allow for commenting out experiments
            if l.startswith('#') or not l:
                continue
            submitted_experiments.add(l.strip())

    print ('There are ' + str(len(submitted_experiments)) + ' experiments')

    # phase 2 - in one pass, obtain lists of object accessions to submit
    output_lists = {'biosamples': [], 'controls': [], 'files': [], 
                    'donors': {'HumanDonor': [], 'FlyDonor': [], 'WormDonor': [], 'MouseDonor': []}, 
                    'experiments': []}
 
    # For purposes of obtaining accessions
    rep_fields = ['library.biosample.accession', 'library.biosample.status', 
                       'library.biosample.donor', 'library.biosample.donor.status',
                       'library.biosample.donor.@type', 'status']

    donor_embedded_fields = ['organism.scientific_name', 'lab.title', 'references',
                                  'documents', 'characterizations']
                                  
    characterization_fields = ['@id', 'characterization_method', 'references.identifiers']

    document_fields = ['document_type', 'attachment.href',
                            '@id', 'urls', 'references']

    gm_fields = ['@id', 'purpose', 'category', 'introduced_tags', 'reagents',
                      'method', 'zygosity', 'modified_site_by_coordinates',
                      'modified_site_by_target_id.label', 
                      'modified_site_by_target_id.investigated_as', 'treatments']

    treatment_fields = ['@id', 'treatment_type','dbxrefs','treatment_term_name','treatment_term_id',
                             'amount','amount_units','duration','duration_units', 'temperature',
                             'temperature_units']

    # top level fields and species-specific fields
    values_to_retain = {
    	'Homo sapiens': ['url','accession','sex','age','age_units','life_stage','health_status',
                         'ethnicity','dbxrefs', 'siblings', 'children', 'twin', 'twin_type'],
    	'Mus musculus': ['url','accession','strain_name','strain_background','dbxrefs', 'source', 
                         'genetic_modifications'],
    	'Drosophila melanogaster': ['url','accession','strain_name','strain_background', 'genotype',
                                    'dbxrefs', 'source', 'genetic_modifications'],
    	'Caenorhabditis elegans': ['num_times_outcrossed','outcrossed_strain.accession','url',
                                   'accession','strain_name','strain_background','genotype',
                                   'dbxrefs', 'source', 'genetic_modifications']}

    # has format dataframe_column_name: output_object_name
    donor_fields_direct_conversion = {'Organism': 'organism', # organsim.scientific_name
                               'Lab': 'lab', # lab.title
                               'URL': 'url',
                               'Accession': 'accession', 'Strain name': 'strain_name',
                               'Strain background': 'strain_background',	
                               'Strain genotype': 'genotype',
                               'Sex': 'sex', 'Age units': 'age_units',
                               'Life stage': 'life_stage', 'Health status': 'health_status',
                               'Ethnicity': 'ethnicity', 
                               'Twin': 'twin', 'Twin type': 'twin_type',
                               'Number of times outcrossed': 'num_times_outcrossed',
                               'outcrossed_strain.accession': 'outcrossed_strain'}

    # fields that need preprocessing 
    donor_fields_indirect_conversion = {'External resources': 'dbxrefs',
                                        'References': 'references',
                                        'Source': 'source', 'Documents': 'documents', 
                                        'Characterizations': 'characterizations', 
                                        'Genetic modifications': 'genetic_modifications',
                                        'Children': 'children', 'Siblings': 'siblings',
                                        'Age': 'age',}

    document_fields_conversion = {'ID': '@id', 'Document type': 'document_type', 
                                 'attachment.href': 'attachment', 'URLs': 'urls', 
                                 'References': 'references'}

    characterizations_fields_conversion = {'ID': '@id', 'Method': 'characterization_method',
                                           'references.identifiers': 'references'}

    gm_fields_conversion = {'ID': '@id', 'Purpose': 'purpose', 'Category': 'category', 
                            'Introduced protein tags': 'introduced_tags', 'Reagents': 'reagents',
                            'Method': 'method', 'Modification zygosity': 'zygosity', 
                            'Modification site': 'modified_site_by_coordinates',
                            'modified_site_by_target_id.label': 'modified_site_by_target_id', 
                            'Treatments': 'treatments'}
                                        
    donor_fields_conversion = {'direct': donor_fields_direct_conversion, 
                                    'indirect': donor_fields_indirect_conversion,
                                    'documents': document_fields_conversion, 
                                    'donor_characterization': characterizations_fields_conversion,
                                    'genetic_modification': gm_fields_conversion}                     

    donor_fields = {}
    donor_fields['HumanDonor'] = donor_embedded_fields + values_to_retain['Homo sapiens']
    donor_fields['MouseDonor'] = donor_embedded_fields + values_to_retain['Mus musculus'] 
    donor_fields['FlyDonor'] = donor_embedded_fields + values_to_retain['Drosophila melanogaster']
    donor_fields['WormDonor'] = donor_embedded_fields + values_to_retain['Caenorhabditis elegans']

    donor_fields = {'documents': document_fields, 'genetic_modification': gm_fields, 
                   'donor_characterization': characterization_fields, 'donor': donor_fields,
                   'treatment': treatment_fields}

    # Get all biosample, donor, control, and file accessions
    get_all_accessions(submitted_experiments, output_lists, rep_fields)

    # Organizes subobjects by type to facilitate retrieval
    donor_objs_by_type = {'document': {}, 'genetic_modification': {}, 'donor_characterization': {}, 
                          'treatment': {}, 'publication': {}}
    
    # Holds all of the subobjects, indexed by @id
    donor_objs_by_id = {}
    
    donor_objs = build_donors(output_lists, donor_fields, donor_fields_conversion, donor_objs_by_type, 
                                     donor_objs_by_id)

    biosample_fields = ['summary', 'accession', 'biosample_ontology.classification', 'biosample_ontology.term_name', 
                        'biosample_ontology.term_id', 'description', 'dbxrefs', 'passage_number',
                        'model_organism_mating_status', 'subcellular_fraction_term_name',
                        'subcellular_fraction_term_id', 'phase', 'url', 'fly_synchronization_stage',
                        'post_synchronization_time', 'post_synchronization_time_units',
                        'worm_synchronization_stage', 'age_units', 'sex', 'health_status', 'age', 
                        'life_stage', 'organism.scientific_name', 'lab.title', 'references',
                        'documents', 'characterizations', 'applied_modifications', 'source.title',
                        'donor', 'treatments']
                                  
    biosample_fields_dict = {'documents': document_fields, 'genetic_modification': gm_fields, 
                   'biosample_characterization': characterization_fields, 'biosample': biosample_fields,
                   'treatment': treatment_fields}

    biosample_fields_indirect_conversion = {#'Database external identifiers': 'dbxrefs', 
                                            'External resources': 'dbxrefs',# need to convert to list
                                            'References': 'references', 
                                            'Source': 'source', 'Documents': 'documents', 
                                            'Characterizations': 'characterizations', 
                                            'Applied modifications': 'applied_modifications', 
                                            'Donor': 'donor', 'Treatments': 'treatments',
                                            'Passage number': 'passage_number'}

    biosample_fields_direct_conversion = {'Summary': 'summary', 'Accession': 'accession',
                                          'Biosample classification': 'biosample_type', 'Biosample term name': 'biosample_term_name',
                                          'biosample_ontology.term_id': 'biosample_term_id', 'Description': 'description',
                                          'Model organism mating status': 'model_organism_mating_status',
                                          'Subcellular fraction': 'subcellular_fraction_term_name',
                                          'subcellular fraction term ID': 'subcellular_fraction_term_id',
                                          'Cell cycle phase': 'phase', 'URL': 'url',
                                          'Fly synchronization stage': 'fly_synchronization_stage',
                                          'Post-synchronization time': 'post_synchronization_time',
                                          'Post-synchronization time units': 'post_synchronization_time_units',
                                          'Worm synchronization stage': 'worm_synchronization_stage',
                                          'Age units': 'age_units', 'Sex': 'sex', 'Health status': 'health_status',
                                          'Age': 'age', 'Life stage': 'life_stage', 'Lab': 'lab',
                                          'Source': 'source', 'Organism': 'organism'}
                                           
    biosample_fields_conversion = {'direct': biosample_fields_direct_conversion, 
                                        'indirect': biosample_fields_indirect_conversion, 
                                        'documents': document_fields_conversion, 
                                        'biosample_characterization': characterizations_fields_conversion,
                                        'genetic_modification': gm_fields_conversion} 

    biosample_objs_by_type = {'document': {}, 'genetic_modification': {}, 'biosample_characterization': {}, 
                              'treatment': {}, 'publication': {}}
    # Holds all of the subobjects, indexed by @id
    biosample_objs_by_id = {}
    biosample_objs = build_biosamples(output_lists, biosample_fields_dict, biosample_fields_conversion, 
                                             biosample_objs_by_type, biosample_objs_by_id)

    experiment_fields = ['date_released', 'accession', 'biosample_ontology.classification', 'assay_title', 'assay_term_name', 'assembly', 
                         'description', 'dbxrefs', 'biosample_ontology.term_name', 'replicates', 'files', 
                         'lab.title', 'references', 'documents', 'possible_controls.accession', 
                         'target.investigated_as', 'target.label']
    
    file_fields = ['alternate_accessions', 'status', 'paired_end', 'derived_from', 'paired_with', 'assembly', 
                   'genome_annotation', 'accession', 'md5sum', 'output_type','file_format', 'file_type', 'href',
                   'content_md5sum', 'read_length', 'read_length_units', 'file_size', 'run_type', 'output_category',
                   'replicate.biological_replicate_number', 'replicate.technical_replicate_number', 'platform.dbxrefs',
                   'platform.term_name', '@id', 'platform', 'replicate', 'index_of']

    excluded_file_output_types = ['filtered reads']

    replicate_fields = ['biological_replicate_number','technical_replicate_number', 'status', '@id']

    library_fields = ['accession','nucleic_acid_starting_quantity_units','nucleic_acid_term_name', 'extraction_method', 
                      'fragmentation_methods', 'library_size_selection_method','size_range','nucleic_acid_starting_quantity', 
                      'spikeins_used', 'documents', 'biosample.accession', 'status', 'biosample.status', '@id']

    replicate_fields.extend('library.' + field for field in library_fields)
    spikein_fields = ['@id', 'accession','dbxrefs', 'description']
    experiment_fields_dict = {'documents': document_fields, 'experiment': experiment_fields,
                              'file': file_fields, 'replicate': replicate_fields, 'reference': spikein_fields}
                         
    experiment_fields_indirect_conversion = {'External identifiers': 'dbxrefs', 
                                             'References': 'references', 'Documents': 'documents', 
                                             'Replicates': 'replicates', 'Files': 'files', 
                                             'possible_controls.accession': 'possible_controls',
                                             'Target of assay': 'target', 'Genome assembly': 'assembly'}
                                             
    experiment_fields_direct_conversion = {'Date released': 'date_released', 'Accession': 'accession', 
                                           'Biosample classification': 'biosample_type', 'Assay title': 'assay_title',
                                           'Assay name': 'assay_term_name',  'Description': 'description', 
                                           'Biosample term name': 'biosample_term_name', 'Lab': 'lab'}
                                           
    file_fields_conversion = {'Alternate accessions': 'alternate_accessions', 'Status': 'status', 
                              'Paired end identifier': 'paired_end', 'Derived from': 'derived_from', 
                              'Paired with': 'paired_with', 'Genome assembly': 'assembly', 
                              'Genome annotation': 'genome_annotation', 'Accession': 'accession', 
                              'MD5sum': 'md5sum', 'Output type': 'output_type', 'File Format': 'file_format', 
                              'File type': 'file_type', 'Download URL': 'href', 'Content MD5sum': 'content_md5sum',
                              'Read length': 'read_length', 'Read length units': 'read_length_units', 
                              'File size': 'file_size', 'Run type': 'run_type', 
                              'Output category': 'output_category', 
                              'replicate.biological_replicate_number': 'biological_replicate_number', 
                              'replicate.technical_replicate_number': 'technical_replicate_number',
                              'platform.dbxrefs': 'dbxrefs', 'platform.term_name': 'term_name',
                              'ID': '@id', 'Platform': 'platform', 'Replicate': 'replicate', 'Index of': 'index_of'}
                              
    replicate_fields_conversion = {'Biological replicate': 'biological_replicate_number', 
                                   'Technical replicate': 'technical_replicate_number', 'ID': '@id',
                                   'Status': 'status'}

    # Add in the library fields to conversion dict
    replicate_fields_conversion.update(zip(('library.' + field for field in library_fields), library_fields))

    experiment_fields_conversion = {'direct': experiment_fields_direct_conversion, 
                                         'indirect': experiment_fields_indirect_conversion, 
                                         'documents': document_fields_conversion, 
                                         'file': file_fields_conversion,
                                         'replicate': replicate_fields_conversion}

    experiment_objs_by_type = {'file': {}, 'replicate': {}, 'document': {}, 'publication': {}, 'reference': {}}

    # Holds all of the subobjects, indexed by @id
    experiment_objs_by_id = {}
    experiment_objs = build_experiments(output_lists, experiment_fields_dict, experiment_fields_conversion, 
                                               experiment_objs_by_type, experiment_objs_by_id)

    exclude_undesired_file_output_types(excluded_file_output_types, experiment_objs, output_lists) 

    # Write to jsons
    write_to_json('donor', donor_objs)
    write_to_json('biosample', biosample_objs)
    write_to_json('experiment', experiment_objs)

    # Obtain s3 file paths
    write_files_file(output_lists['files'], splitext(args.infile)[0] + '_FILES_TO_UPLOAD.txt')

    t1 = time.process_time()
    print('elapsed time (seconds): ', t1 - t0)


def encoded_get(url, keypair=None, frame='object', return_response=False):
    """
    This is used to obtain s3 file paths by write_files_file
    """
    url_obj = urlsplit(url)
    new_url_list = list(url_obj)
    query = parse_qs(url_obj.query)
    if 'format' not in query:
        new_url_list[3] += '&format=json'
    if 'frame' not in query:
        new_url_list[3] += '&frame=%s' % (frame)
    if 'limit' not in query:
        new_url_list[3] += '&limit=all'
    if new_url_list[3].startswith('&'):
        new_url_list[3] = new_url_list[3].replace('&', '', 1)
    get_url = urlunsplit(new_url_list)
    max_retries = 10
    max_sleep = 10
    while max_retries:
        try:
            if keypair:
                response = requests.get(get_url, auth=keypair, headers=GET_HEADERS)
            else:
                response = requests.get(get_url, headers=GET_HEADERS)
        except (requests.exceptions.ConnectionError,
                requests.exceptions.SSLError) as e:
            print >> sys.stderr, e
            sleep(max_sleep - max_retries)
            max_retries -= 1
            continue
        else:
            if return_response:
                return response
            else:
                return response.json()


def getKeyPair(path_to_key_pair_file, server_name):
    keysf = open(path_to_key_pair_file, 'r')
    keys_json_string = keysf.read()
    keysf.close()
    keys = json.loads(keys_json_string)
    key_dict = keys[server_name]
    AUTHID = key_dict['key']
    AUTHPW = key_dict['secret']
    return (AUTHID, AUTHPW)


def write_to_json(object_type, objects):
    """
    Write experiment, biosample, and donor objects to json to different 
    subdirectories of a folder. Log accessions to submission_script.log
    """
    object_type = object_type + 's'
    print('STARTING {}'.format(object_type).upper())
    accessions = []
    for obj in objects:
        accession = obj['accession']
        print(accession)
        accessions.append(accession)
        # file_path = "../{}/".format(object_type) + accession + "_modified.json"
        # file_path = "/Users/paul/geo-debug/new_script/{}/".format(object_type) + accession + "_modified.json"
        file_path = "../{}/{}_modified.json".format(object_type, accession)
        with open(file_path, "w") as file_out:
            file_out.write(json.dumps(obj, indent=4, sort_keys=True))
    
    # Log accessions
    with open('submission_script.log', 'a') as f:
        f.write(strftime('%Y-%m-%d %H:%M'))
        f.write(' Outputted {}\n'.format(object_type))
        for accession in accessions:
            f.write(accession+'\n')
        
    print('FINISHED {}'.format(object_type).upper())


def get_all_accessions(submitted_experiments, output_lists, rep_fields):
    """
    Obtains control, biosample, donor, and file accessions from submitted_experiments.
    The addition of these values to output_lists depends on the status of the 
    objects on the portal (generally, if they are released).
    """
    print('Obtaining accessions of objects to be submitted')
    get_experiment_accessions(submitted_experiments, output_lists, rep_fields)             
    get_experiment_accessions(set(output_lists['controls']), output_lists, rep_fields)
    for donor_type, donors in output_lists['donors'].items():
        output_lists['donors'][donor_type] = list(set(donors))


def make_list_from_str(list_as_str):
    """ 
    Handles report fields coming from json arrays, including NaNs and comma-
    separated values.
    """
    if not list_as_str or pd.isna(list_as_str):
        return []
    else:
        return list_as_str.strip('[]').split(',')


def get_experiment_accessions(accession_set, output_lists, rep_fields):
    """
    Get all experiments, then get a report on the replicates, extracting controls, biosamples, 
    and donors accessions to output_lists.
    """
    if not accession_set:
        return
    released_experiments = []
    replicates_by_exp = {}
    response = search('Experiment', list(accession_set))
    for exp in response:
        if exp['status'] in ['released', 'archived']:
            # Store experiment object for later use
            output_lists['experiments'].append(exp)
            released_experiments.append(exp['accession'])
            # extract controls, checking if not empty
            if exp['possible_controls']:
                control_accessions = [control_id.split('/')[2] for control_id in exp['possible_controls']]
                # output_lists['controls'].extend(exp['possible_controls'])
                output_lists['controls'].extend(control_accessions)
            # extract replicates
            replicates = exp['replicates']
            replicates_by_exp[exp['accession']] = replicates
    
    all_replicates = list(set(v for sublist in replicates_by_exp.values() for v in sublist))
    # get biosamples and donors from all replicates in one pass
    rep_df = get_report_tsv_from_fields('replicate', all_replicates, rep_fields)
    for j in rep_df.index:
        if rep_df.loc[j, 'library.biosample.status'] in ['released', 'archived']:
            biosample_accessions = rep_df.loc[j, 'library.biosample.accession']
            output_lists['biosamples'].extend(make_list_from_str(biosample_accessions))
            if rep_df.loc[j, 'library.biosample.donor.status'] in ['released', 'archived']:
                donors = make_list_from_str(rep_df.loc[j, 'library.biosample.donor'])
                donor_type = rep_df.loc[j, 'library.biosample.donor.@type'].split(',')[0]
                for donor in donors:
                    donor_accession = donor.split('/')[2]
                    output_lists['donors'][donor_type].append(donor_accession)
    # get fastq file accessions from the released experiments
    get_file_accessions(list(set(released_experiments)), output_lists)


def get_file_accessions(experiments, output_lists):
    # Search for files in the experiment
    for experiment in experiments:
        res = search('File', experiment)
        # Check for criteria
        for file in res:
            if file['status'] in ('uploading', 'released', 'in progress', 'archived') and not_in_sra(file):
                if 'restricted' in file:
                    assert file['restricted'] == False, 'Aborting, cannot submit restricted file {}'.format(file['accession'])
                output_lists['files'].append(file['accession'])
                    # output_lists['files'].append(file['s3_uri'])


def build_experiments(output_lists, fields_dict, fields_conversion, objs_by_type, objs_by_id):
    print('Building experiment objects to submit')
    experiment_fields_list = fields_dict['experiment']
    
    # holds experiment objects to be submitted
    objs_to_output = []

    # first pass: collect fields that we get ordered information from report
    experiments = [exp['accession'] for exp in output_lists['experiments']]
    if not experiments:
        return objs_to_output

    experiment_df = get_report_tsv_from_fields('Experiment', experiments, experiment_fields_list)

    # collect ids of properties we need to get separately
    for i in experiment_df.index:
        experiment_obj = minimize_experiment(experiment_df, i, fields_conversion)
        objs_to_output.append(experiment_obj)
        objs_by_type['document'][experiment_df.loc[i, 'Accession']] = make_list_from_str(experiment_df.loc[i, 'Documents'])
        objs_by_type['file'][experiment_df.loc[i, 'Accession']] = make_list_from_str(experiment_df.loc[i, 'Files'])
        objs_by_type['replicate'][experiment_df.loc[i, 'Accession']] = make_list_from_str(experiment_df.loc[i, 'Replicates'])
        objs_by_type['publication'][experiment_df.loc[i, 'Accession']] = make_list_from_str(experiment_df.loc[i, 'References'])

    # replicate.library has documents that are aggregated with experiment documents
    get_objects_w_report('replicate', objs_by_type, fields_dict, fields_conversion, objs_by_id)
    
    get_objects('document', objs_by_type, objs_by_id)

    # Obtain publication @ids from documents 
    for documents in objs_by_type['document'].values():
        for document in documents:
            objs_by_type['publication'][document] = objs_by_id[document]['references']
    
    all_document_ids = [value for vals in objs_by_type['document'].values() for value in vals]
    all_documents = [objs_by_id[doc_id] for doc_id in all_document_ids]
    
    get_objects('publication', objs_by_type, objs_by_id)
    
    # Backfill the documents with the publications (into the reference field)
    # backfill_fields(objs_by_type['document'].values(), objs_by_id, is_document=True)
    backfill_fields(all_documents, objs_by_id, is_document=True)
    
    get_objects('reference', objs_by_type, objs_by_id)
    
    # Backfill spikeins_used in libraries
    backfill_fields(objs_by_type['replicate'].values(), objs_by_id)
    
    get_objects_w_report('file', objs_by_type, fields_dict, fields_conversion, objs_by_id)
    
    # Add replicate objects to objs_by_id for backfilling into experiments
    for k, v in objs_by_type['replicate'].items():
        objs_by_id[k] = v 
    
    # Backfill the files, replicates, and documents
    backfill_fields(objs_to_output, objs_by_id, is_experiment=True)
    
    return objs_to_output


def minimize_experiment(report_df, report_df_index, fields_conversion):
    """
    Takes an experiment dataframe and minimizes the row given by report_df_index.
    The main feature of this function compared to other minimizers is the 
    unique handling of the targets.
    """
    obj = {}
    fields_direct_conversion = fields_conversion['direct']
    fields_indirect_conversion = fields_conversion['indirect']
    for column in report_df.columns:
        # columns whose values can be directly submitted
        if column in fields_direct_conversion:
            # don't submit fields with blank (NaN) values except for description
            if pd.isna(report_df.loc[report_df_index, column]):
                if column == 'Description':
                    obj[fields_direct_conversion[column]] = ''
                continue
            obj[fields_direct_conversion[column]] = report_df.loc[report_df_index, column]
        # columns that need preprocessing
        elif column in fields_indirect_conversion:
            if column == 'Target label':
                if pd.isna(report_df.loc[report_df_index, column]):
                    pass
                else:
                    # Process target depending on if invesigated as control
                    investigated_as = make_list_from_str(report_df.loc[report_df_index, 'target.investigated_as'])
                    if 'control' in investigated_as:
                        obj[fields_indirect_conversion[column]] = 'Control'
                    else:
                        obj[fields_indirect_conversion[column]] = report_df.loc[report_df_index, column]
            else:
                l = make_list_from_str(report_df.loc[report_df_index, column])
                obj[fields_indirect_conversion[column]] = l
    
    return obj


def build_biosamples(output_lists, fields_dict, fields_conversion, objs_by_type, objs_by_id):
    """
    Takes a list of biosample accessions and gives list of biosample objects (dicts) for submission.
    Unlike build_donors, there are no subtypes of biosample, so this function
    is a little simpler than build_donors.
    """
    print('Building biosample objects to submit')
    biosample_fields_list = fields_dict['biosample']
    
    # holds donor objects to be submitted
    objs_to_output = []

    # first pass: collect fields that we get ordered information from report
    biosamples = output_lists['biosamples']
    if not biosamples:
        return objs_to_output

    biosample_df = get_report_tsv_from_fields('Biosample', biosamples, biosample_fields_list)

    # collect ids of properties we need to get separately
    for i in biosample_df.index:
        biosample_obj = minimize_donors_and_biosamples(biosample_df, i, fields_conversion)
        objs_to_output.append(biosample_obj)
        objs_by_type['document'][biosample_df.loc[i, 'Accession']] = make_list_from_str(biosample_df.loc[i, 'Documents'])
        objs_by_type['genetic_modification'][biosample_df.loc[i, 'Accession']] = make_list_from_str(biosample_df.loc[i, 'Applied modifications'])
        objs_by_type['biosample_characterization'][biosample_df.loc[i, 'Accession']] = make_list_from_str(biosample_df.loc[i, 'Characterizations'])
        objs_by_type['treatment'][biosample_df.loc[i, 'Accession']] = make_list_from_str(biosample_df.loc[i, 'Treatments'])
        objs_by_type['publication'][biosample_df.loc[i, 'Accession']] = make_list_from_str(biosample_df.loc[i, 'References'])

    # get info that report has trouble handling
    get_objects('document', objs_by_type, objs_by_id)
    get_objects('publication', objs_by_type, objs_by_id)
    # need references.identifiers from donor_characterization report
    get_objects_w_report('biosample_characterization', objs_by_type, fields_dict, fields_conversion, objs_by_id)
    # need things like modifications_by_site_id.investigated_as 
    get_objects_w_report('genetic_modification', objs_by_type, fields_dict, fields_conversion, objs_by_id)
    # get treatment objects corresponding to gms
    # for some reason these are unicode strings?
    get_objects('treatment', objs_by_type, objs_by_id)
    # backfill gm objects with treatments 
    backfill_fields(objs_by_type['genetic_modification'].values(), objs_by_id)
    # add backfilled gms to objs_by_id
    for k, v in objs_by_type['genetic_modification'].items():
        objs_by_id[k] = v 
    # backfill all donor objects with gms, documents, and characterizations
    backfill_fields(objs_to_output, objs_by_id)

    return objs_to_output


def build_donors(output_lists, report_fields, fields_conversion, objs_by_type, objs_by_id):
    """
    Takes a list of donors and gives list of donor objects (dicts) for submission.
    """
    print('Building donor objects to submit')
    donor_fields = report_fields['donor']
    
    # holds donor objects to be submitted
    objs_to_output = []

    desired_species = ('Homo sapiens', 'Mus musculus', 'Drosophila melanogaster', 'Caenorhabditis elegans')
    # first pass: collect fields that we get ordered information from report 
    for donor_type, donors in output_lists['donors'].items():
        if not donors:
            continue
        donor_fields_list = donor_fields[donor_type]
        donor_df = get_report_tsv_from_fields(donor_type, donors, donor_fields_list)
        # donor_df.to_csv('/Users/paul/donor_fields_debug.tsv', sep = '\t')
        # collect ids of properties we need to get separately
        for i in donor_df.index:
            donor_accession = donor_df.loc[i, 'Accession']
            if donor_df.loc[i, 'Organism'] not in desired_species:
                continue
            else:
                donor_obj = minimize_donors_and_biosamples(donor_df, i, fields_conversion)
                objs_to_output.append(donor_obj)
            objs_by_type['document'][donor_accession] = make_list_from_str(donor_df.loc[i, 'Documents'])
            if donor_type != 'HumanDonor':
                objs_by_type['genetic_modification'][donor_accession] = make_list_from_str(donor_df.loc[i, 'Genetic modifications'])
            objs_by_type['donor_characterization'][donor_accession] = make_list_from_str(donor_df.loc[i, 'Characterizations'])
            objs_by_type['publication'][donor_accession] = make_list_from_str(donor_df.loc[i, 'References'])

    # get info that report has trouble handling
    get_objects('document', objs_by_type, objs_by_id)
    get_objects('publication', objs_by_type, objs_by_id)
    # need references.identifiers from donor_characterization report
    get_objects_w_report('donor_characterization', objs_by_type, report_fields, fields_conversion, objs_by_id)
    # need things like modifications_by_site_id.investigated_as 
    get_objects_w_report('genetic_modification', objs_by_type, report_fields, fields_conversion, objs_by_id)
    # get treatment objects corresponding to gms
    get_objects('treatment', objs_by_type, objs_by_id)
    # backfill gm objects with treatments 
    backfill_fields(objs_by_type['genetic_modification'].values(), objs_by_id)
    # add backfilled gms to objs_by_id
    for k, v in objs_by_type['genetic_modification'].items():
        objs_by_id[k] = v 
    # backfill all donor objects with gms, documents, and characterizations
    backfill_fields(objs_to_output, objs_by_id)

    return objs_to_output


def minimize_donors_and_biosamples(report_df, report_df_index, fields_conversion):
    """
    This function generically converts the report tsvs into submittable objects,
    at least for biosamples and donors. Empty (NaN) fields are generally not included 
    in the output objects. However, documents, references, characterizations, and 
    genetic_modifications (for nonhuman donors) always appear in the output object as 
    a list even if they are empty.
    """
    # initialize output object
    obj_to_output = {}
    fields_direct_conversion = fields_conversion['direct']
    fields_indirect_conversion = fields_conversion['indirect']
    for column in report_df.columns:
        value = report_df.loc[report_df_index, column]
        # columns whose values can be directly submitted
        if column in fields_direct_conversion:
            # don't submit fields with blank (NaN) values
            if pd.isna(value):
                if column == 'Description':
                    obj_to_output[fields_direct_conversion[column]] = ''
                continue
            elif column == 'Age':
                obj_to_output[fields_direct_conversion[column]] = str(value)
            else:
                obj_to_output[fields_direct_conversion[column]] = value
        # columns that need preprocessing
        elif column in fields_indirect_conversion:
            if column in ('Source', 'Donor'):
                # check if source is NaN
                if pd.isna(value):
                    continue
                # convert from '/sources/source-name/' to 'source-name'
                parsed = value.split('/')[2]
                obj_to_output[fields_indirect_conversion[column]] = parsed
            elif column == 'Age':
                obj_to_output[fields_indirect_conversion[column]] = str(value)
            elif column == 'Passage number':
                if not pd.isna(value):
                    obj_to_output[fields_indirect_conversion[column]] = int(value)
            else:
                l = make_list_from_str(value)
                obj_to_output[fields_indirect_conversion[column]] = l
    
    return obj_to_output 


def backfill_fields(objs_to_output, objs_by_id, is_experiment=False, is_document=False):
    """
    Fields of objects in objs_to_output containing lists of @ids are populated
    with the desired subobjects corresponding to the @ids with the appropriate
    minimizations, if needed. 
    """
    for obj in objs_to_output:
        for field in obj:
            if field == 'documents':
                if is_experiment:
                    obj[field] = [minimize_experiment_document(objs_by_id[obj_id]) for obj_id in obj[field]]
                if not is_experiment:
                    obj[field] = [minimize_document(objs_by_id[obj_id]) for obj_id in obj[field]]
            if field == 'treatments':
                obj[field] = [minimize_treatment(objs_by_id[obj_id]) for obj_id in obj[field]]
            # Need to get to library in replicate
            if field == 'library':
                if 'spikeins_used' in obj[field]:
                    obj[field]['spikeins_used'] = [minimize_reference(objs_by_id[obj_id]) for obj_id in obj[field]['spikeins_used']]
                if 'documents' in obj[field]: 
                    obj[field]['documents'] = [minimize_experiment_document(objs_by_id[obj_id]) for obj_id in obj[field]['documents']]
            if field == 'references' and is_document:
                # We need a dictionary of references in format @id:[identifier(s)]
                obj[field] = {obj_id:minimize_publication(objs_by_id[obj_id]) for obj_id in obj[field]}
            if field == 'references' and not is_document:
                obj[field] = [minimize_publication(objs_by_id[obj_id]) for obj_id in obj[field]]
            if field == 'replicates':
                # handle situation where replicate is not released and not in objs_by_id
                output = []
                for obj_id in obj[field]:
                    if obj_id not in objs_by_id:
                        continue
                    else:
                        output.append(objs_by_id[obj_id])
                obj[field] = output
            if field in ('characterizations', 'genetic_modifications', 'applied_modifications', 
                         'files'):
                # these objects have already been minimized
                obj[field] = [objs_by_id[obj_id] for obj_id in obj[field]]


def minimize_publication(publication):
    return publication['identifiers']


def minimize_reference(reference):
    """
    Takes a reference spikein dataset object and strips it down to a subset of its fields.
    """
    reference_fields = ['accession', 'dbxrefs', 'description']
    minimized_reference = {field:reference[field] for field in reference_fields if field in reference}
    return minimized_reference


def minimize_experiment_document(document):
    """
    Takes a document belonging to an experiment or to a library in an experiment
    and strips it down to a subset of desired fields. This differs from other 
    non-experiment documents in that the attachment is a dictionary rather than
    a simple string concatenation of document @id and hrefself.
    """
    minimized_document = {}
    for key in ('document_type', 'urls', 'references', 'attachment'):
        if key in document:
            if key == 'attachment':
                minimized_document[key] = minimize_attachment(document[key], document['@id'])
            else:
                minimized_document[key] = document[key]

    return minimized_document


def minimize_attachment(attachment, doc_id):
    """
    Takes an attachment obtained from a document in an experiment object and 
    strips it down to a subset of desired fields. The document @id, given by
    doc_id, is prepended to the attachment href.
    """
    minimized_attachment = {}
    
    for key in ('md5sum', 'href'):
        if key in attachment:
            if key == 'href':
                minimized_attachment[key]=doc_id + attachment[key]
            elif key == 'md5sum':
                minimized_attachment[key] = attachment[key]
    
    return minimized_attachment


def minimize_document(document):
    """
    Takes a document obtained directly from its json from the portal and strips
    it down to a subset of desired fields. The document @id is prepended to the
    attachment href
    """
    minimized_document = {}
    
    for field in ('document_type', 'urls', 'references'):
        if field in document:
            minimized_document[field] = document[field]
    
    if 'attachment' in document:
        minimized_document['attachment'] = document['@id'] + document['attachment']['href']

    return minimized_document


def minimize_treatment(treatment):
    """
    Takes a treatment obtained directly from its json from the portal and strips
    it down to a subset of desired fields.
    """
    minimized_treatment = {}
    
    for key in ['treatment_type','dbxrefs','treatment_term_name','treatment_term_id','concentration',
                'concentration_units','duration','duration_units','temperature','temperature_units']:
        if key in treatment.keys():
            minimized_treatment[key]=treatment[key] 
    
    return minimized_treatment


def minimize_gms(gm_df, fields_conversion, objs_by_type, objs_by_id):
    """
    Takes a dataframe containing gms and prepares objects containing most gm properties 
    (except treatments, which are backfilled later). Both gms and treatments are 
    added to objs_by_id, since they are both backfilled into other objects later. 
    """
    # gm_list = objs_by_type['genetic_modification']
    gm_fields_conversion = fields_conversion['genetic_modification']
    
    for i in gm_df.index:
        minimized_gm = {}
        for column in gm_df.columns:
            value = gm_df.loc[i, column]
            if pd.isna(value):
                # Filter out fields with no value
                continue
            if column in ['Purpose', 'Category', 'Method','Modification zygosity']:
                # Insert field/value pair into output object that don't need conversion
                minimized_gm[gm_fields_conversion[column]] = value
            if column == 'Modification site':
                # The raw value is given by the report as a string corresponding to a json object
                processed = ast.literal_eval(value)
                minimized_gm[gm_fields_conversion[column]] = processed
            if column in ('Introduced protein tags', 'Reagents'):
                # The raw data is a string corresponding to a list of json objects,
                # although without enclosing brackets
                raw_data = '[' + value + ']'
                processed = ast.literal_eval(raw_data)
                minimized_gm[gm_fields_conversion[column]] = processed
            if column == 'Treatments':
                # Needs conversion to a list
                values_list = make_list_from_str(value)
                gm_accession = gm_df.loc[i, 'ID'].split('/')[2]
                # add treatments to objects to get
                objs_by_type['treatment'][gm_accession] = values_list
                minimized_gm[gm_fields_conversion[column]] = values_list
            if column == 'modified_site_by_target_id.label':
                # Check if control target, if it is not, add target label to gm object
                if gm_df.loc[i, 'modified_site_by_target_id.investigated_as'] != 'control':
                    minimized_gm[gm_fields_conversion[column]] = value
        # objs_by_id[gm_df.loc[i, 'ID']] = minimized_gm
        objs_by_type['genetic_modification'][gm_df.loc[i, 'ID']] = minimized_gm


def minimize_files(file_df, fields_conversion, objs_by_type, objs_by_id):
    """This is for the purpose of minimizing files to incorporate into experiment
    objects, not to submit to GEO"""
    file_fields_conversion = fields_conversion['file']
    for i in file_df.index:
        file_dict = {}
        platform = {}
        replicate =  {}
        for column in file_df.columns:
            value = file_df.loc[i, column]
            if column == 'Alternate accessions':
                # This always appears in output file object, need to check first before filtering out NaNs 
                values_list = make_list_from_str(value)
                file_dict[file_fields_conversion[column]] = values_list
            elif column == 'platform.dbxrefs':
                # Want to add this to platform to avoid being filtered out if NaN
                platform[file_fields_conversion[column]] = make_list_from_str(value)
            elif pd.isna(value) or column in ('ID', 'Platform', 'Replicate'):
                continue
            elif column in ('File size', 'Read length'):#, 'Paired End Identifier'):
                file_dict[file_fields_conversion[column]] = int(value)
            elif column == 'Paired end identifier':
                # In original script, value is outputted as string
                file_dict[file_fields_conversion[column]] = str(value).rstrip('.0')
            elif column == 'Derived from':
                values_list = make_list_from_str(value)
                file_dict[file_fields_conversion[column]] = values_list
            elif column == 'Index of':
                values_list = make_list_from_str(value)
                file_dict[file_fields_conversion[column]] = values_list
            elif column == 'platform.term_name':
                platform[file_fields_conversion[column]] = value
            elif column in ('replicate.biological_replicate_number', 'replicate.technical_replicate_number'):
                replicate[file_fields_conversion[column]] = int(value)
            else:
                file_dict[file_fields_conversion[column]] = value
            # elif column == 'platform.dbxrefs':
            #     platform[file_fields_conversion[column]] = make_list_from_str(file_df.loc[i, column])

        if not pd.isna(file_df.loc[i, 'Platform']):
            file_dict['platform'] = platform
        if not pd.isna(file_df.loc[i, 'Replicate']):
            file_dict['replicate'] = replicate

        objs_by_id[file_df.loc[i, 'ID']] = file_dict
        
        
def minimize_replicates(rep_df, fields_conversion, objs_by_type, objs_by_id):
    """
    Convert dataframe comtaining rows with replicate and associated library 
    metadata into a replicate object containing a library subobject.
    """
    rep_fields_conversion = fields_conversion['replicate']
    for i in rep_df.index:
        if rep_df.loc[i, 'Status'] not in ['released', 'archived']:
            continue
        replicate = {}
        library = {}
        library_accession = rep_df.loc[i, 'library.@id'].split('/')[2]
        for column in rep_df.columns:
            value = rep_df.loc[i, column]
            if column == 'library.spikeins_used':
                # Always want spikeins list in library, even if empty (NaN in report)
                values = make_list_from_str(value)
                library[rep_fields_conversion[column]] = values
                objs_by_type['reference'][library_accession] = values
            elif column == 'library.documents':
                values = make_list_from_str(value)
                library[rep_fields_conversion[column]] = values
                objs_by_type['document'][library_accession] = values
            elif pd.isna(value):# or column in ('ID', 'Platform', 'Replicate'):
                continue
            elif column in ('Biological replicate', 'Technical replicate'):
                replicate[rep_fields_conversion[column]] = int(value)
            elif column.startswith('library'):
                if column in ('library.biosample.status', 'library.@id', 'library.status'):
                    continue
                elif column == 'library.biosample.accession':
                    if rep_df.loc[i, 'library.biosample.status'] in ['released', 'archived']:
                        library['biosample'] = value
                elif column == 'library.nucleic_acid_starting_quantity':
                    library[rep_fields_conversion[column]] = str(int(value))
                elif column == 'library.fragmentation_methods':
                    fragmentation_methods = make_list_from_str(value)
                    library[rep_fields_conversion[column]] = ', '.join(fragmentation_methods)
                else:
                    library[rep_fields_conversion[column]] = value
            else:
                # Other columns don't make it into outputted experiment json
                pass
        
        if rep_df.loc[i, 'library.status'] in ['released', 'archived']:
            replicate['library'] = library

        objs_by_type['replicate'][rep_df.loc[i, 'ID']] = replicate


def get_objects(object_type, objs_by_type, objs_by_id):
    """
    Given a list of @ids for desired objects (obtained from objs_by_type using the object_type key), 
    obtain the objects via search using frame=object and output them to a dictionary.
    """
    # for documents and characterizations
    documents_dict = objs_by_type[object_type]
    all_documents = [doc for docs in documents_dict.values() for doc in docs]
    all_documents = list(set(all_documents))
    if all_documents:
        res = search(object_type, all_documents)
        for obj in res:
            objs_by_id[obj['@id']] = obj


def get_objects_w_report(object_type, objs_by_type, fields_dict, fields_conversion, objs_by_id):
    """
    Given a list of @ids for desired objects (obtained from objs_by_type using the object_type key), 
    obtain the desired fields for all objects using a report tsv and process them into objects, 
    outputting them into the objs_by_id dictionary where the object's @id is the key
    """
    # for documents and characterizations
    documents_dict = objs_by_type[object_type]
    document_fields = fields_dict[object_type]
    document_fields_conversion = fields_conversion[object_type]
    all_documents = [doc for docs in documents_dict.values() for doc in docs]
    all_documents = list(set(all_documents))

    if all_documents:
        output_columns = list(document_fields_conversion.keys())
        output_columns.remove('ID')
        # output_columns = document_fields_conversion.keys().remove('ID')
        doc_df = get_report_tsv_from_fields(object_type, all_documents, document_fields)
        # can't put directly into objs_by_id since treatment is not filled in
        if object_type == 'genetic_modification':
            minimize_gms(doc_df, fields_conversion, objs_by_type, objs_by_id)
            return
        if object_type == 'file':
            minimize_files(doc_df, fields_conversion, objs_by_type, objs_by_id)
            return
        if object_type == 'replicate':
            minimize_replicates(doc_df, fields_conversion, objs_by_type, objs_by_id)
            return
        for i in doc_df.index:
            # Essentially create a minimized output object, really only for characterizations
            document_obj = {}
            for column in output_columns:
                if column in doc_df.columns:
                    # These fields need to be treated as lists
                    if column in ('URLs', 'References', 'references.identifiers'):
                        values_list = make_list_from_str(doc_df.loc[i, column])
                        document_obj[document_fields_conversion[column]] = values_list
                    else:
                        document_obj[document_fields_conversion[column]] = doc_df.loc[i, column]
            objs_by_id[doc_df.loc[i, 'ID']] = document_obj


def capitalize_object_type(object_type):
    split = object_type.split("_")
    if len(split) == 1:
        if "Donor" in object_type:
            return object_type
        return object_type.capitalize()
    return "".join(i.capitalize() for i in split)


def search(object_type, accession):
    # Max nginx url
    max_url_length = 8192
    object_type = capitalize_object_type(object_type)
    url_part_0 = 'search/?type='
    url_part_2 = '&frame=object&limit=all'
    dataset_key = '&dataset='
    id_key = '&@id='
    accession_key = '&accession='
    file_format_key = '&file_format=fastq'
    
    fixed_url_components = SERVER + url_part_0 + object_type + url_part_2
    fixed_length = len(fixed_url_components)
    available_url = max_url_length - fixed_length
    
    # search for fastqs in a given dataset, no need to chunk since always short query
    if object_type == 'File':
        # search_key = dataset_key
        # # searching for files by dataset requires extra parameter
        # available_url -= len(file_format_key)
        experiment_id = '/experiments/{}/'.format(accession)
        url = SERVER + url_part_0 + object_type + dataset_key + experiment_id + file_format_key + url_part_2
        response = requests.get(url, auth=(AUTHID, AUTHPW))
        res = response.json()
        return res['@graph']
        
    elif object_type in ('Experiment', 'Biosample'):
        search_key = accession_key
    else:
        search_key = id_key
    
    # In case of searching for @id, need to account for replacement of '@' symbol
    search_key_long = search_key.replace('@', '%40')
    # Compute number of items allowed per search, use integer division to floor
    chars_per_item = len(search_key_long + accession[0].replace('/', '%2f').replace('@', '%40'))
    # Take floor via integer division to ensure below url size limit
    num_items = available_url // chars_per_item
    # Chunk @id list according to num_items to allow search
    chunked = [accession[i:i + num_items] for i in range(0, len(accession), num_items)]
    
    responses = []
    
    for chunk in chunked:
        url_part_1 = url_part_0 + object_type + search_key + search_key.join(chunk)
        URL = SERVER + url_part_1 + url_part_2
        response = requests.get(URL, auth=(AUTHID, AUTHPW))
        res = response.json()
        responses.extend(res['@graph'])
    return responses
    

def get_report_tsv_from_fields(object_type, accession, fields_list):
    """
    Given an accession or a list of accessions, return a dataframe containing
    the report results. The max URL allowed by nginx is 8192 characters. 
    / and @ characters are converted to %2f and %40, respectively, which 
    greatly increases the url size. Therefore, reports are generate in chunks
    then concatenated.
    """
    # Max nginx url
    max_url_length = 8192
    object_type = capitalize_object_type(object_type)
    url_part_0 = 'report.tsv?type='
    id_key = '&@id='
    accession_key = '&accession='
    # This is automatically appended to report queries, but need to check url size
    limit_param = '&limit=all'
    
    # Need to know how long this is to chunk the reports
    url_part_2 = build_report_fields_url(fields_list)
    # Needed for @id
    url_part_2_long = url_part_2.replace('@', '%40')
    
    fixed_url_components = SERVER + url_part_0 + object_type + url_part_2_long + limit_param
    fixed_length = len(fixed_url_components)
    available_url = max_url_length - fixed_length
    
    if isinstance(accession, list):
        # Figure out whether searching @ids or accessions
        if object_type in ('Replicate', 'GeneticModification', 'DonorCharacterization',
                           'BiosampleCharacterization', 'File'):
            search_key = id_key
        else: 
            search_key = accession_key
        
        # In case of searching for @id, need to account for replacement of '@' symbol
        search_key_long = search_key.replace('@', '%40')
        # Compute number of items allowed per search, use integer division to floor
        chars_per_item = len(search_key_long + accession[0].replace('/', '%2f').replace('@', '%40')) 
        num_items = available_url // chars_per_item
        # Chunk @id list according to num_items to allow search
        chunked = [accession[i:i + num_items] for i in range(0, len(accession), num_items)]
        dfs = []
        # Obtain dfs corresponding to the report on each chunk of accessions
        for chunk in chunked:
            url_part_1 = url_part_0 + object_type + search_key + search_key.join(chunk)
            URL = SERVER + url_part_1 + url_part_2
            response = requests.get(URL, auth=(AUTHID, AUTHPW))
            res_df = pd.read_csv(pd.compat.StringIO(response.text), sep = '\t', header=1)
            dfs.append(res_df)
        # Combine all of the df, ignoring index to avoid conflicting indices
        all_dfs = pd.concat(dfs, ignore_index=True, sort=True)
        return all_dfs
    else:
        url_part_1 = url_part_0 + object_type + accession_key + accession

    URL = SERVER + url_part_1 + url_part_2
    response = requests.get(URL, auth=(AUTHID, AUTHPW))
    res_df = pd.read_csv(pd.compat.StringIO(response.text), sep = '\t', header=1)
    return res_df


def build_report_fields_url(fields_list):
    specified_fields = ['&field=' + field for field in fields_list]
    return ''.join(specified_fields) + '&limit=all'


def not_in_sra(file_object):
    if 'dbxrefs' not in file_object or len(file_object['dbxrefs']) == 0:
        return True
    else:
        sra_flag = False
        for entry in file_object['dbxrefs']:
            if entry.startswith('SRA:'):
                return False
        if not sra_flag:
            return True


def write_files_file(files_to_upload, file_path):
    print('STARTING FILES')
    with open(file_path, 'w') as file_of_files:
        # for s3_uri in set(files_to_upload):
        for file_accession in set(files_to_upload):
            file_object = encoded_get(SERVER+'/files/'+file_accession, keypair)
            s3_path_url = file_object.get("s3_uri")
            file_of_files.write('/s3'+s3_path_url[4:]+'\n')
    print('FINISHED FILES')


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'infile',
        help='input file containing accessions of experiments to be submitted.'
    )
    args = parser.parse_args()
    return args

def exclude_undesired_file_output_types(excluded_file_output_types, experiment_objs, output_lists):
    file_accessions_being_removed = set()
    for current_experiment in experiment_objs:
        # set up temporary list of files to keep
        files_to_keep = []
        for current_file in current_experiment['files']:
            current_file_output_type = current_file['output_type']
            if current_file_output_type in excluded_file_output_types:
                current_file_accession = current_file['accession']
                file_accessions_being_removed.add(current_file_accession)
            else:
                files_to_keep.append(current_file)
        # update the files to keep with the experiment
        current_experiment['files'] = files_to_keep   
    # loop through files and exclude files being removed
    new_file_list = []
    old_file_set = set(output_lists['files'])
    for current_file_accession in old_file_set:
        if current_file_accession not in file_accessions_being_removed:
            new_file_list.append(current_file_accession)
    output_lists['files'] = new_file_list

# Globals
HEADERS = {'accept': 'application/json'}
GET_HEADERS = {'accept': 'application/json'}
POST_HEADERS = {'accept': 'application/json',
                'Content-Type': 'application/json'}
#SERVER = "https://test.encodedcc.org/"
SERVER = 'https://www.encodeproject.org/'

keypair = getKeyPair('keypairs.json', 'test')

AUTHID = keypair[0]
AUTHPW = keypair[1]

if __name__ == '__main__':
    main()
