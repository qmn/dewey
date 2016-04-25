Dewey
=====
Dewey is intended to be a *fast* placer and router for Minecraft circuits. It
is the C port of [PERSHING](https://github.com/qmn/pershing). An external logic
synthesis program, like [yosys](http://www.clifford.at/yosys/), is needed to
synthesize Verilog into a Berkeley Logic Interchange File (BLIF). `dewey` will
take such a file and produce a compacted layout of logic cells and redstone
wires needed to produce a functioning circuit in Minecraft.

Requirements
------------
- `make`
- a C compiler (e.g., `gcc`)

Why is it called Dewey?
--------------------------
PERSHING is named after a [ballistic missle
system](https://en.wikipedia.org/wiki/MGM-31_Pershing), which is, in turn,
named after an American general, John Pershing. I wanted to find Pershing's
peer in the U.S. Navy, on which I settled on George Dewey (the "C" version of
Pershing).

Other Notes
-----------
A publication resulting from the creation of PERSHING, which is the original
version of this project, appeared at the first annual conference on TBD
([SIGTBD](http://sigtbd.csail.mit.edu/)), a joke conference at MIT.
