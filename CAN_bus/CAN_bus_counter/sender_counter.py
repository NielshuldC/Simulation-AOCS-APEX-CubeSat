#main.py -- put your code here!S
import pyb 
from pyb import CAN
led = pyb.LED(3)

button = pyb.Pin(pyb.Pin.board.SW, pyb.Pin.IN)

def debounce(pin):
    cur_value = pin.value()
    active = 0
    while active < 20:
        if pin.value() != cur_value:
            active += 1
        else:
            active = 0
        pyb.delay(1)  

while True:
    debounce(button)
    can = CAN(1, CAN.NORMAL, extframe=True, prescaler=16, sjw=4, bs1=25, bs2=1)
    can.setfilter(0, CAN.LIST16, 0, (123, 124, 125, 126))
    can.send('1', 123)
    pyb.delay(100)
    led.toggle()
    pyb.delay(100)
