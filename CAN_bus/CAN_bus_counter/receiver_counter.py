#main.py -- put your code here!
import pyb 
from pyb import CAN
led = pyb.LED(3)


pinbit = [pyb.Pin(pyb.Pin.board.D49, pyb.Pin.OUT),pyb.Pin(pyb.Pin.board.D50, pyb.Pin.OUT),pyb.Pin(pyb.Pin.board.D65, pyb.Pin.OUT),pyb.Pin(pyb.Pin.board.D64, pyb.Pin.OUT)]

intval = 0

for i in range(0, 4):
    pinbit[i].value(1)
pyb.delay(500)
for i in range(0, 4):
    pinbit[i].value(0)
pyb.delay(500)

while True:
    try:
        can = CAN(1, CAN.NORMAL, extframe=True, prescaler=16, sjw=4, bs1=25, bs2=1)
        can.setfilter(0, CAN.LIST16, 0, (123, 124, 125, 126))
        (idr, isRTR, filterMatchIndex, telegram) = can.recv(0) 
        if idr == 123:
            binval = bin(intval)
            lenght = len(binval) - 2
            state = [0,0,0,0]
            if lenght == 1:
                state[0] = 0
                state[1] = 0
                state[2] = 0
                state[3] = binval[2]
                for i in range(0, 4):
                    pinbit[i].value(int(state[i]))
                pyb.delay(100)
                intval = intval + 1
                pyb.delay(100)
            elif lenght == 2:
                state[0] = 0
                state[1] = 0
                state[2] = binval[2]
                state[3] = binval[3]
                for i in range(0, 4):
                    pinbit[i].value(int(state[i]))
                pyb.delay(100)
                intval = intval + 1
                pyb.delay(100)
            elif lenght == 3:
                state[0] = 0
                state[1] = binval[2]
                state[2] = binval[3]
                state[3] = binval[4]
                for i in range(0, 4):
                    pinbit[i].value(int(state[i]))
                pyb.delay(100)
                intval = intval + 1
                pyb.delay(100)
            elif lenght == 4:
                state[0] = binval[2]
                state[1] = binval[3]
                state[2] = binval[4]
                state[3] = binval[5]
                for i in range(0, 4):
                    pinbit[i].value(int(state[i]))
                pyb.delay(100)
                intval = intval + 1
                pyb.delay(100)
            else:
                intval=0
                state = [0,0,0,0]
                for i in range(0, 4):
                    pinbit[i].value(int(state[i]))
                pyb.delay(100)
                intval = intval + 1
                pyb.delay(100)   
        else:
            print("NO addition to the counter/")       
    except:
        led.toggle()
