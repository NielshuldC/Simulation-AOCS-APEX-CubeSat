import pyb 
from pyb import CAN
led = pyb.LED(3)
serial=pyb.USB_VCP()
lines=[]
nbLines=4

while(True):
    try:
        can = CAN(1, CAN.NORMAL, extframe=True, prescaler=16, sjw=4, bs1=25, bs2=1)
        can.setfilter(0, CAN.LIST16, 0, (123, 122, 125, 126))
        (idr, isRTR, filterMatchIndex, telegram) = can.recv(0)
        if idr == 122:
            print(telegram)
       

    while (serial.any() and len(lines)<=nbLines) :
        line=serial.readline()
        if "br=" in line or "tr=" in line:
            line = line[:-2]
            data = line.split("=")[1]
            line = str(data)
        lines.append(line) 
    if len(lines)==nbLines :
        can = CAN(1, CAN.NORMAL, extframe=True, prescaler=16, sjw=4, bs1=25, bs2=1)
        can.setfilter(0, CAN.LIST16, 0, (123, 122, 125, 126))
        can.send(lines, 123)
        led.toggle()
        lines=[]
