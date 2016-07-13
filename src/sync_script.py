'''
Using existing donors, biosamples and experiments files we should create a list of objects to be updated.
Step 1. Using experimental accessions - get all the biosamples and donors and files that supposed to be up in GEO.

Step 2. Run the update using this new list - creating the changes on the branch.

Step 3. Report FASTQ files that have no SRA identifier - suggesting  that they have not been submitted yet to GEO. 
It is not reliable - because files may ave SRA identifiers for different reason. The only way I see is to keep list of files we submitted already and check against it?
Or use the github capabilities to detect changes - and go through updated experiments looking for any additional fastqs - and creating a list of these.
Probably that is cleaner. The assumption beneath is that any file up untill now is submitted. (what is not updated - was submitted.)

'''
import glob
import requests
import json


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
HEADERS = {'accept': 'application/json'}


# phase 1 - collect all experiments submitted so far.

submittedExperiments = set()
for filename in glob.glob('../experiments/*.json'):
    submittedExperiments.add(filename.split('/')[2].split('_')[0])

# phase 2 - go over the experiments submitted so far and create a set of biosamples and donors 
controls_list = []
biosamples_list = []
mone = 0
for experiment in submittedExperiments:
    mone += 1
    URL = "https://www.encodeproject.org/" + experiment + \
          "/?frame=embedded&format=json"
    response = requests.get(URL, auth=(AUTHID, AUTHPW), headers=HEADERS)
    experiment = response.json()
    print (str(mone) + ' Inspecting Experiment ')  # + str(experiment))
    controls_list.extend(extract_controls(experiment))
    biosamples_list.extend(extract_biosamples(experiment))
    if mone > 5:
        break

biosample_objects = []
mone = 0
for biosample in biosamples_list:
    mone += 1
    accession_number = biosample
    URL = "https://www.encodeproject.org/"+accession_number+"/?frame=embedded&format=json"
    response = requests.get(URL, auth=(AUTHID, AUTHPW), headers=HEADERS)
    biosample_object = response.json()
    print str(mone) + ' Inspecting Biosample ' + str(accession_number)
    biosample_objects.append(biosample_object)


donors_list = extract_donors(biosample_objects)

print (set(controls_list))
print (set(biosamples_list))
print (set(donors_list))

# phase 3 - go over the experiments submitted and get the list of FASTQs


files_list = []
mone = 0
for experiment in submittedExperiments:
    mone += 1
    URL = "https://www.encodeproject.org/search/?type=File&dataset=/experiments/" + experiment + \
          "/&file_format=fastq&format=json&frame=object"
    response = requests.get(URL, auth=(AUTHID, AUTHPW), headers=HEADERS)
    all_experiment_fastqs = response.json()
    experimental_fastqs = [f for f in all_experiment_fastqs['@graph'] if f['status'] not in ['deleted', 'revoked', 'replaced', 'upload failed', 'format check failed', 'archived']]
    for fastq_file in experimental_fastqs:
        acc = fastq_file['accession']
        files_list.append(acc)
    if mone > 5:
        break
print (files_list)