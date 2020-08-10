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
            state = [0,0,0,0]
            if telegram >= 3.5 and telegram <= 4.5 :
                state[0] = 0
                state[1] = 1
                state[2] = 1
                state[3] = 0
                for i in range(0, 4):
                    pinbit[i].value(int(state[i]))
                pyb.delay(100)
            elif telegram < 3.5:
                state[0] = 1
                state[1] = 1
                state[2] = 0
                state[3] = 0
                for i in range(0, 4):
                    pinbit[i].value(int(state[i]))
                pyb.delay(100)
            elif telegram > 4.5:
                state[0] = 0
                state[1] = 0
                state[2] = 1
                state[3] = 1
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
            print("NO position data/")       
    except:
        led.toggle()
