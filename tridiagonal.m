   f   u   n   c   t   i   o   n           [   V   ,   V   i   n   v   ]       =       t   r   i   d   i   a   g   o   n   a   l   (   n   )      
                   a       =       r   a   n   d   n   (   n   ,   1   )   ;      
                   b       =       r   a   n   d   n   (   n   -   1   ,   1   )   ;      
                   c       =           r   a   n   d   n   (   n   -   1   ,   1   )   ;      
                      
                   V       =       (   d   i   a   g   (   a   )   +   d   i   a   g   (   b   ,   1   )   +   d   i   a   g   (   c   ,   -   1   )   )   ;      
                   T   h   e   t   a       =       z   e   r   o   s   (   n   +   1   ,   1   )   ;      
                   T   h   e   t   a   (   1   )       =       1   ;       T   h   e   t   a   (   2   )       =       a   (   1   )   ;      
                   P   h   i       =       z   e   r   o   s   (   n   +   1   ,   1   )   ;      
                   P   h   i   (   n   +   1   )       =   1   ;       P   h   i   (   n   )       =       a   (   n   )   ;      
                      
                   V   i   n   v       =       z   e   r   o   s   (   n   )   ;      
                   f   o   r       i       =       2   :   n      
                                   T   h   e   t   a   (   i   +   1   )       =       a   (   i   )   *   T   h   e   t   a   (   i   )       -       b   (   i   -   1   )   *   c   (   i   -   1   )   *   T   h   e   t   a   (   i   -   1   )   ;      
                                   j       =       n   +   1   -   i   ;      
                                   P   h   i   (   j   )       =       a   (   j   )   *   P   h   i   (   j   +   1   )       -       b   (   j   )   *   c   (   j   )   *   P   h   i   (   j   +   2   )   ;      
                   e   n   d      
                      
                   f   o   r       j       =       1   :   n      
                                   f   o   r       i       =       1   :   j   -   1      
                                                   V   i   n   v   (   i   ,   j   )       =       (   -   1   )   ^   (   i   +   j   )   *       p   r   o   d   (   b   (   i   :   j   -   1   )   )   *       T   h   e   t   a   (   i   )   *   P   h   i   (   j   +   1   )   /   T   h   e   t   a   (   n   +   1   )   ;      
                                   e   n   d      
                                   V   i   n   v   (   j   ,   j   )       =       T   h   e   t   a   (   j   )   *   P   h   i   (   j   +   1   )   /   T   h   e   t   a   (   n   +   1   )   ;      
                                   f   o   r       i       =       j   +   1   :   n      
                                                   V   i   n   v   (   i   ,   j   )       =       (   -   1   )   ^   (   i   +   j   )   *       p   r   o   d   (   c   (   j   :   i   -   1   )   )   *       T   h   e   t   a   (   j   )   *   P   h   i   (   i   +   1   )   /   T   h   e   t   a   (   n   +   1   )   ;      
                                   e   n   d      
                   e   n   d