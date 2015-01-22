## jellyfish-mirror ##
For lack of a better name, this tool will use a binary .jf file (as outputted by Jellyfish) to then count another set of input files.
It will count only k-mers which appear in the given .jf file and uses the same hash function, so `jellyfish dump` will output them in the same sorted order.


## todo ##
 - does not support auto dumping and merging of files, as jellyfish-count does

## notes ##
this is not feature equivalent to `jellyfish count` for things like read filtering etc.
