import glob
import requests
import json
import ExperimentBoiler
import geoDonorMinimiser
import geoBiosampleMinimiser
import urlparse
import sys
from time import sleep


HEADERS = {'accept': 'application/json'}
GET_HEADERS = {'accept': 'application/json'}
POST_HEADERS = {'accept': 'application/json',
                'Content-Type': 'application/json'}
#SERVER = "https://test.encodedcc.org/"
SERVER = "https://www.encodeproject.org/"

def encoded_get(url, keypair=None, frame='object', return_response=False):
    url_obj = urlparse.urlsplit(url)
    new_url_list = list(url_obj)
    query = urlparse.parse_qs(url_obj.query)
    if 'format' not in query:
        new_url_list[3] += "&format=json"
    if 'frame' not in query:
        new_url_list[3] += "&frame=%s" % (frame)
    if 'limit' not in query:
        new_url_list[3] += "&limit=all"
    if new_url_list[3].startswith('&'):
        new_url_list[3] = new_url_list[3].replace('&', '', 1)
    get_url = urlparse.urlunsplit(new_url_list)
    max_retries = 10
    max_sleep = 10
    while max_retries:
        try:
            if keypair:
                response = requests.get(get_url,
                                        auth=keypair,
                                        headers=GET_HEADERS)
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


def extract_biosamples(exp):
    samples = []
    if exp['status'] == 'released' and \
       'replicates' in exp and \
       len(exp['replicates']) > 0:
        for replicate in exp['replicates']:
            if replicate['status'] == 'released' and \
               replicate['library']['status'] == 'released' and \
               replicate['library']['biosample']['status'] == 'released':
                samples.append(replicate['library']['biosample']['accession'])
    return list(set(samples))


def extract_controls(exp):
    if "possible_controls" in exp and \
       len(exp['possible_controls']) > 0:
        controls_list = []
        for e in exp['possible_controls']:
            controls_list.append(e['accession'])

        return list(set(controls_list))
    else:
        return []


def extract_donors(biosamples_list):
    donors = []
    for biosample in biosamples_list:
        if biosample['status'] == 'released' and \
           'donor' in biosample and \
           biosample['donor']['status'] == 'released':
            donors.append(biosample['donor']['accession'])
    return list(set(donors))


keypair = getKeyPair('keypairs.json', 'test')

AUTHID = keypair[0]
AUTHPW = keypair[1]

# phase 1 - collect all experiments submitted so far.

submittedExperiments = set()

exp_f = open('ggr.10.12.2017', 'r')

#exp_f = open("try_exp_list", "r")
for l in exp_f:
    submittedExperiments.add(l.strip())
exp_f.close()

print ('There are ' + str(len(submittedExperiments)) + ' experiments')

# phase 2 - go over the experiments submitted so far and create a set of biosamples and donors 
controls_list = []
biosamples_list = []
mone = 0
for experiment in submittedExperiments:
    mone += 1
    URL = SERVER + experiment + "/?frame=embedded&format=json"
    response = requests.get(URL, auth=(AUTHID, AUTHPW), headers=HEADERS)
    experiment_o = response.json()
    print (str(mone) + ' Inspecting Experiment ')  # + str(experiment))
    if experiment_o['status'] == 'released':
        controls_list.extend(extract_controls(experiment_o))
        biosamples_list.extend(extract_biosamples(experiment_o))
    #if mone > 5:
    #    break

for experiment in set(controls_list):
    mone += 1
    URL = SERVER + experiment + "/?frame=embedded&format=json"
    response = requests.get(URL, auth=(AUTHID, AUTHPW), headers=HEADERS)
    experiment_o = response.json()
    print (str(mone) + ' Inspecting Control Experiment ')  # + str(experiment))
    if experiment_o['status'] == 'released':
        biosamples_list.extend(extract_biosamples(experiment_o))
    #if mone > 5:
    #    break


biosample_objects = []
mone = 0
for biosample in biosamples_list:
    mone += 1
    accession_number = biosample
    #print ('biosample accession' + accession_number)
    URL = SERVER+accession_number+"/?frame=embedded&format=json"
    response = requests.get(URL, auth=(AUTHID, AUTHPW), headers=HEADERS)
    biosample_object = response.json()
    print str(mone) + ' Inspecting Biosample ' + str(accession_number)
    biosample_objects.append(biosample_object)


donors_list = extract_donors(biosample_objects)


# phase 3 - go over the experiments submitted and get the list of FASTQs


files_list = []
mone = 0

experiments_and_controls = submittedExperiments | set(controls_list)
released_experiments = set()
for experiment in experiments_and_controls:
    mone += 1
    print str(mone) + ' Inspecting Experiment ' + str(experiment)
    URL = SERVER + experiment + "/?frame=embedded&format=json"
    response = requests.get(URL, auth=(AUTHID, AUTHPW), headers=HEADERS)
    experiment_o = response.json()
    if experiment_o['status'] == 'released':
        released_experiments.add(experiment_o['accession'])
        mone += 1
        URL = SERVER + "search/?type=File&dataset=/experiments/" + experiment + \
            "/&file_format=fastq&format=json&frame=object"
        response = requests.get(URL, auth=(AUTHID, AUTHPW), headers=HEADERS)
        all_experiment_fastqs = response.json()
        experimental_fastqs = [f for f in all_experiment_fastqs[
            '@graph'] if f['status'] not in [
            'deleted', 'revoked', 'replaced',
            'upload failed', 'format check failed', 'content error',
            'archived']]
        for fastq_file in experimental_fastqs:
            acc = fastq_file['accession']
            files_list.append(acc)
        #if mone > 5:
        #    break

files_to_upload = []
mone = 0
for f in files_list:
    mone += 1
    URL = SERVER + f + "/?frame=embedded&format=json"
    response = requests.get(URL, auth=(AUTHID, AUTHPW), headers=HEADERS)
    file_object = response.json()
    if 'dbxrefs' not in file_object or len(file_object['dbxrefs']) == 0:
        files_to_upload.append(f)
    else:
        sra_flag = False
        for entry in file_object['dbxrefs']:
            if entry.startswith('SRA:') is True:
                sra_flag = True
        if sra_flag is False:
            files_to_upload.append(f)

#print (set(files_to_upload))
# print (set(experiments_and_controls))
# print (set(biosamples_list))
#print (set(donors_list))
 
'''
 at this point we have the following:
 (1) donors_list
 (2) biosamples_list
 (3) submittedExperiments that contains the controls as well
 (3) files_to_upload - a list of file accessions that need to be uploaded to SRA
'''


for donor_accession in set(donors_list):
    print (donor_accession)
    URL = SERVER+donor_accession+"/?frame=embedded&format=json"
    response = requests.get(URL, auth=(AUTHID, AUTHPW), headers=HEADERS)
    response_json_dict = response.json()
    file_out = open("../donors/" + donor_accession+"_modified.json", "w")
    file_out.write((json.dumps(geoDonorMinimiser.minimise_donor(response_json_dict), indent=4, sort_keys=True)))
    file_out.close()
print ('FINISHED DONORS')

for biosample_accession in set(biosamples_list):
    print (biosample_accession)
    URL = SERVER+biosample_accession+"/?frame=embedded&format=json"
    response = requests.get(URL, auth=(AUTHID, AUTHPW), headers=HEADERS)
    response_json_dict = response.json()
    file_out = open("../biosamples/" + biosample_accession+"_modified.json", "w")
    file_out.write((json.dumps(geoBiosampleMinimiser.minimise_biosample(response_json_dict), indent=4, sort_keys=True)))
    file_out.close()
print ('FINISHED BIOSAMPLES')


print ('STARTING EXPERIMENTS')
for experimental_accession in released_experiments:
    print (experimental_accession)
    URL = SERVER+experimental_accession+"/?frame=embedded&format=json"
    response = requests.get(URL, auth=(AUTHID, AUTHPW), headers=HEADERS)
    response_json_dict = response.json()
    file_out = open("../experiments/" + experimental_accession+"_modified.json", "w")
    file_out.write((json.dumps(ExperimentBoiler.boildown_experiment(response_json_dict), indent=4, sort_keys=True)))
    file_out.close()

print ('FINISHED EXPERIMENTS')


print ('STARTING FILES')
file_of_files = open('NEW_FILES_OCTOBER_2017_GGR_TO_UPLOAD', 'w')



for file_accession in set(files_to_upload):
    up_creds = encoded_get(SERVER+'/files/'+file_accession+'/@@upload', keypair)
    s3_path_url = up_creds['@graph'][0]['upload_credentials']['upload_url']
    file_of_files.write('/s3'+s3_path_url[4:]+'\n')
file_of_files.close()

print ('FINISHED FILES')