import pyb
import time
import logging; logging.get_logger('spacecan').set_level(logging.DEBUG)
import spacecan

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


MY_NODE_ID = 0  # master node has id 0
SLAVE_NODE_IDS = [123]  # list of all attached slave nodes

bus_a, bus_b = spacecan.Bus(1), spacecan.Bus(2)
network = spacecan.Network(MY_NODE_ID, bus_a, bus_b)
network.set_filters([  # receive telemetry frames from all slaves
    (spacecan.ID_TM + node_id, spacecan.FULL_MASK)
    for node_id in SLAVE_NODE_IDS])
heartbeat = spacecan.HeartbeatProducer(network)
sync = spacecan.SyncProducer(network)
telecommand = spacecan.TelecommandProvider(network)
telemetry = spacecan.TelemetryConsumer(network)

network.start()
heartbeat.start(500)
sync.start(5000)


def main():
    while True:
        debounce(button)
            print("Send telecommand, mode:", mode)
            telecommand.send(SLAVE_NODE_IDS[0], mode)
            mode = "OFF" if mode == "ON" else "ON"

main()  # run the main loop

