
# Python3 program to merge overlapping Intervals
# in O(n Log n) time and O(1) extra space
def mergeIntervals(arr):

        # Sorting based on the increasing order
        # of the start intervals
        arr.sort(key = lambda x: x[0])

        # array to hold the merged intervals
        m = []
        s = -10000
        max = -100000
        for i in range(len(arr)):
            a = arr[i]
            if a[0] > max:
                if i != 0:
                    m.append([s,max])
                max = a[1]
                s = a[0]
            else:
                if a[1] >= max:
                    max = a[1]

        #'max' value gives the last point of
        # that particular interval
        # 's' gives the starting point of that interval
        # 'm' array contains the list of all merged intervals

        if max != -100000 and [s, max] not in m:
            m.append([s, max])

        return m
        # print("The Merged Intervals are :", end = " ")
        # for i in range(len(m)):
        #     print(m[i], end = " ")

# main
if __name__ == '__main__':
    arr = [[6, 8], [1, 9], [2, 4], [4, 7]]
    print(mergeIntervals(arr))
