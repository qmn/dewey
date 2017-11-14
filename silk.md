input:
- dimensions of design
- matrix of (currently) occupied spaces (i.e., extracted cell placements)
- locations of pins of a given net, and their direction

procedure:
1. create a min-heap for each pin in the net
2. create a union-find structure for the pins
3. initialize `n_groups` to len(pins)
4. initialize a distance matrix, doing the following:
 1. mark each pin (or pin + partial net) with `BT_START` and distance 0

5. initialize a visited matrix, doing the following:
 1. mark all locations as unvisited
 2. mark all current groups

6. while `n_groups` > 1:
 1. select the group with the lowest min-heap element
 2. expand the heap as before
 3. if the lowest min-heap element expands into a group not currently in my union-find set:
  1. create a new "segment" with these two points as endpoints
  2. with the intersecting point, backtrace to the remote `BT_START`
  3. with the point from which we expanded, backtrace to the local `BT_START`
  4. add this segment with the addition of the two backtraces
  5. delete the min-heap corresponding to both sets
  3. mark all points along the partial net as visited, and add these points to the new min-heap for expansion

OLD
  1. backtrace to the nearest `BT_START`, generating a new partial net
  2. mark all parts of the to-be-joined groups as unvisited, marking each pin and partial net with `BT_START` and distance 0
  3. join the two groups in the union-find structure
  4. decrement `n_groups`

output:
- list of blocks connecting all pins of net
