f = open('files_trans.3.14.17', 'r')
m = open('files.3.14.17', 'w')
for l in f:
    # print ('line: ' + l)
    arr = l.strip().split()
    if len(arr) == 2:
        m.write(arr[1].split('.')[0] + '\t' + arr[0]+'\n')
    elif len(arr) == 3:
        m.write(arr[1].split('.')[0] + '\t' + arr[0]+'\n')
        m.write(arr[2].split('.')[0] + '\t' + arr[0]+'\n')
f.close()
m.close()
