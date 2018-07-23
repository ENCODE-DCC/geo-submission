f = open('files_trans.5.2.18', 'r')
m = open('files.5.2.18', 'w')
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
