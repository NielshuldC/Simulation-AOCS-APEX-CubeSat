import pyb 
from pyb import CAN
led = pyb.LED(3)

serial=pyb.USB_VCP()

while (serial.any() and len(lines)<=nbLines)
    line=serial.readline() 
    can = CAN(1, CAN.NORMAL, extframe=True, prescaler=16, sjw=4, bs1=25, bs2=1)
    can.setfilter(0, CAN.LIST16, 0, (123, 124, 125, 126))
    can.send(line, 123)
    pyb.delay(100)
    led.toggle()
    pyb.delay(100)
