#main.py -- put your code here!
# Sending message board
import pyb 
from pyb import CAN
led = pyb.LED(3)# Using red LED from board to check if messages are sent
while True:
    led.toggle()
    pyb.delay(500)
    can = CAN(1, CAN.NORMAL, extframe=True, prescaler=16, sjw=4, bs1=25, bs2=1)
    # Receiver and sender should have the same parameters for the CAN bus
    can.setfilter(0, CAN.LIST16, 0, (123, 124, 125, 126)) 
    can.send('1', 123)
    # sending message '1' with ID 123
    pyb.delay(500)


 
