# main.py -- put your code here!
# Message Receiving board
# Code for testing the CAN bus when not screening your board
import pyb
from pyb import CAN # importing CAN libraries
led = pyb.LED(3) # Using red LED on board to monitor when CAN messages are not received
aled = pyb.LED(2) # Using blue LED on board to monitor when CAN messages are received
while True:
    try:
        can = CAN(1, CAN.NORMAL, extframe=True, prescaler=16, sjw=4, bs1=25, bs2=1)
        # extframe=True for 29-Bit CAN ID 
        # Baudrate = Clock_frequency/[prescaler*(1+bs1+bs2)]
        # to check zour frequency of your clock:
        # freq = pyb.freq()
        # freq
        can.setfilter(0, CAN.LIST16, 0, (123, 124, 125, 126))
        # Setting filter is compulsory for receiving messages
        (idr, isRTR, filterMatchIndex, telegram) = can.recv(0) 
        if idr == 123 : #ID of the sender
            aled.toggle()# Blue light toggles - message received
        else:
            print("NO messages received check CAN bus/")       
    except:
        led.toggle() # Red light toggles - no message- error on CAN bus
