# Set the task name
TASK = obsvis_fov

# Uncomment the correct choice indicating either SKA or TST flight environment
FLIGHT_ENV = SKA

BIN = obsvis_fov
SHARE = fov.py

include /proj/sot/ska/include/Makefile.FLIGHT

.PHONY: dist install docs version

install: 
#  Uncomment the lines which apply for this task
	mkdir -p $(INSTALL_BIN)
	mkdir -p $(INSTALL_SHARE)
#	mkdir -p $(INSTALL_DATA)
#	mkdir -p $(INSTALL_DOC)

	rsync --times --cvs-exclude $(BIN) $(INSTALL_BIN)/
	rsync --times --cvs-exclude $(SHARE) $(INSTALL_SHARE)/
#	rsync --times --cvs-exclude $(DATA) $(INSTALL_DATA)/
#	rsync --archive --times --cvs-exclude $(DOC)/ $(INSTALL_DOC)/

