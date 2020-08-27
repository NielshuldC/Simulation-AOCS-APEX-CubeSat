import pyb 
from pyb import CAN
led = pyb.LED(3)

half_period = 2000

def debounce(pin):
    cur_value = pin.value()
    active = 0
    while active < 20:
        if pin.value() != cur_value:
            active += 1
        else:
            active = 0
        pyb.delay(1)  


debounce(button)
print("button pressed")
can = CAN(1, CAN.NORMAL, extframe=True, prescaler=16, sjw=4, bs1=25, bs2=1)
can.setfilter(0, CAN.LIST16, 0, (123, 124, 125, 126))
can.send('1', 122)

i = 0 

while i < 1:
    try:
        can = CAN(1, CAN.NORMAL, extframe=True, prescaler=16, sjw=4, bs1=25, bs2=1)
        can.setfilter(0, CAN.LIST16, 0, (123, 124, 125, 126))
        (idr, isRTR, filterMatchIndex, telegram) = can.recv(0)
        maneuver = telegram[2]
        i = 1

while True: 
        try:
            can = CAN(1, CAN.NORMAL, extframe=True, prescaler=16, sjw=4, bs1=25, bs2=1)
            can.setfilter(0, CAN.LIST16, 0, (123, 124, 125, 126))
            (idr, isRTR, filterMatchIndex, telegram) = can.recv(0)
            if telegram[2] >= maneuver+half_period:
                can = CAN(1, CAN.NORMAL, extframe=True, prescaler=16, sjw=4, bs1=25, bs2=1)
                can.setfilter(0, CAN.LIST16, 0, (123, 124, 125, 126))
                can.send('1', 122)       
    except:
        led.toggle()
