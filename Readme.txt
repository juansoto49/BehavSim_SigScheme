*********************************************************************************************
* This code intends to recreate a behavioral simulation environment                         *
* with different signaling schemes such as NRZ, PAM4, QPSK and QAM16                        *
*                                                                                           *
* This code is free to be modified for needs.                                               *
*                                                                                           *
* Behavioural simulation include TX and RX Equalization and S-Parameter                     *
* pulse respone extraction.                                                                 *
*                                                                                           *
* ------------    ------------    ---------    -------------    -----------    -----------  *
* | Ideal TX |    | Signaling|    | TX FFE|    |S-Parameter|    | RX CTLE |    |Output   |  *
* | Driver   |--->| Scheme   |--->| EQ    |--->| Pulse Resp|--->| + AGC   |--->|Waveform |  *
* |          |    | Mapping  |    |       |    | Extraction|    | + DFE EQ|    |         |  *
* ------------    ------------    ---------    -------------    -----------    -----------  *
*                                                                                           *
* Note: Algorithm doesn't simulate in this order. Simulation order                          *
*                                                                                           *
*                                                                                           *
* 1. Determine channel pulse response                                                       *
* 2. Apply EQ over channel pulse response. Do EQ optimization based on                      *
*	PDA aperture criteria. The flow is recursive and the logic is:                      *
*	a) Apply FFE EQ (set of taps)                                                       *
*		b) Apply CTLE EQ curve (test each curve)                                    *
*			c) Apply DFE to each tap                                            *
* 3. Perform Signaling Scheme Bit Mapping                                                   *
* 4. Create a bit-by-bit behavioral convolution simulation with N random bits               *
*	where the output is a waveform.                                                     *
* 5. Grap the Eye Diagram, PDA, Constellation Diagram and PSD.                              *
*********************************************************************************************