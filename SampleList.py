import sys

FromValue = sys.argv[1]
ToValue = sys.argv[2]
ByValue = sys.argv[3]


for i in range(FromValue, ToValue, ByValue):
    print(i)
