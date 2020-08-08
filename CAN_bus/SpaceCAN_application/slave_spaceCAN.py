import pyb
import logging; logging.get_logger('spacecan').set_level(logging.DEBUG)
import spacecan

led = pyb.LED(3)


pinbit = [pyb.Pin(pyb.Pin.board.D49, pyb.Pin.OUT),pyb.Pin(pyb.Pin.board.D50, pyb.Pin.OUT),pyb.Pin(pyb.Pin.board.D65, pyb.Pin.OUT),pyb.Pin(pyb.Pin.board.D64, pyb.Pin.OUT)]

intval = 0

for i in range(0, 4):
    pinbit[i].value(1)
pyb.delay(500)
for i in range(0, 4):
    pinbit[i].value(0)
pyb.delay(500)

MY_NODE_ID = 123  # must be in range 1 - 127

bus_a, bus_b = spacecan.Bus(1), spacecan.Bus(2)
network = spacecan.Network(MY_NODE_ID, bus_a, bus_b)
network.set_filters([
    # receive sync, hb, and telecommands from master
    (spacecan.ID_SYNC, spacecan.FULL_MASK),
    (spacecan.ID_HEARTBEAT, spacecan.FULL_MASK),
    (spacecan.ID_TC + MY_NODE_ID, spacecan.FULL_MASK)])
heartbeat = spacecan.HeartbeatConsumer(network)
sync = spacecan.SyncConsumer(network)
telecommand = spacecan.TelecommandConsumer(network)
telemetry = spacecan.TelemetryProvider(network)

network.start()
heartbeat.start(500)
sync.start()


def main():
    print("Slave node started")
    mode = "ON"
    while True:
            while telecommand.any():
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


main()  # run the main loop
