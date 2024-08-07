## Charge model training

This directory contains scripts and data involved in fitting a neural network.

"openff-gnn-am1bcc-0.1.0-rc.2.pt" was trained with the following process:

### Pickling training data

Data was initially converted into a [random access format](https://arrow.apache.org/docs/python/ipc.html#writing-and-reading-random-access-files) manually using parallelized workers. This is not strictly
necessary as the nagl library *should* be able to do this on the fly, but doing it this way is firstly
much faster due to the parallel workers, and secondly allows you to control what's going on at a more
fine-grained level.

This was done using `run-prepare-data-pickled.sh`.

### Training the model

