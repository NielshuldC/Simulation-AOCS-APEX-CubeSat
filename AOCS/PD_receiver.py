import pyb 
import math
from pyb import CAN
led = pyb.LED(3)

P = 0.04
D = 0.005
state = 100

while True:
    try:
        can = CAN(1, CAN.NORMAL, extframe=True, prescaler=16, sjw=4, bs1=25, bs2=1)
        can.setfilter(0, CAN.LIST16, 0, (123, 124, 125, 126))
        (idr, isRTR, filterMatchIndex, telegram) = can.recv(0) 
        if idr == 123:

            # Satellite in first quadrant x y
            if telegram[0] != 4 and telegram[1] >= 0 and telegram[2] >= 0 :
            value = (4-telegram[0])*P+(4/state-telegram[3])*D
            if telegram[1] != 0:
                angle = math.atan(telegram[2]/telegram[1])
                line = [value*math.cos(angle),value*math.sin(angle)]
            else:
                line = [0,value]
            can = CAN(1, CAN.NORMAL, extframe=True, prescaler=16, sjw=4, bs1=25, bs2=1)
            can.setfilter(0, CAN.LIST16, 0, (123, 124, 125, 126))
            can.send(line, 122)   

            # Satellite in second quadrant x -y
            elif telegram[0] != 4 and telegram[1] >= 0 and telegram[2] <= 0 :
            value = (4-telegram[0])*P+(4/state-telegram[3])*D
            if telegram[1] != 0:
                angle = math.atan(telegram[2]/telegram[1])
                line = [value*math.cos(angle),-value*math.sin(angle)]
            else:
                line = [0,-value]
            can = CAN(1, CAN.NORMAL, extframe=True, prescaler=16, sjw=4, bs1=25, bs2=1)
            can.setfilter(0, CAN.LIST16, 0, (123, 124, 125, 126))
            can.send(line, 122)   

            # Satellite in third quadrant -x -y
            elif telegram[0] != 4 and telegram[1] <= 0 and telegram[1] <= 0 :
            value = (4-telegram[0])*P+(4/state-telegram[3])*D
            if telegram[1] != 0:
                angle = math.atan(telegram[2]/telegram[1])
                line = [-value*math.cos(angle),-value*math.sin(angle)]
            else:
                line = [0,-value]
            can = CAN(1, CAN.NORMAL, extframe=True, prescaler=16, sjw=4, bs1=25, bs2=1)
            can.setfilter(0, CAN.LIST16, 0, (123, 124, 125, 126))
            can.send(line, 122)  

            # Satellite in forth quadrant -x y
            elif telegram[0] != 4 and telegram[1] <= 0 and telegram[1] >= 0 :
            value = (4-telegram[0])*P+(4/state-telegram[3])*D
            if telegram[1] != 0:
                angle = math.atan(telegram[2]/telegram[1])
                line = [-value*math.cos(angle),value*math.sin(angle)]
            else:
                line = [0,value]
            can = CAN(1, CAN.NORMAL, extframe=True, prescaler=16, sjw=4, bs1=25, bs2=1)
            can.setfilter(0, CAN.LIST16, 0, (123, 124, 125, 126))
            can.send(line, 122)  
            
            else:
                None
         else:
            print("NO position data/")       
    except:
        led.toggle()
