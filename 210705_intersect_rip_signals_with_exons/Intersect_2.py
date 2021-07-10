# 3. Report coordinates of exonic intersection for each peak

file = open("intersect.txt", "r")
lines = file.readlines()
file.close()

new_file = open("intersection_points.txt", "w")
counter = 1
previous = 0
index = 1
for line in lines:
    arr = line.split("\t")
    counter = int(arr[3].strip("Peak_"))
    if counter > index:
        new_file.write('\n')
        index += 1
    if int(arr[7]) == -1:
        new_file.write('\n')
        index += 1
    else:
        if counter > previous:
            new_file.write(str(max(int(arr[1]), int(arr[7]))+1) + "-" + str(min(int(arr[2]), int(arr[8]))))
        elif counter == previous:
            new_file.write(";" + str(max(int(arr[1]), int(arr[7]))+1) + "-" + str(min(int(arr[2]), int(arr[8]))))
    previous = counter
new_file.close()
