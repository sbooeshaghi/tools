To build the project run

```
$ mkdir build; cd build
$ cmake ..
$ make
```

The t2g script can then be ran in one of the following ways:

```
$ cat genes.gtf | t2g make -p - > tr2g.txt
$ t2g make -p - < genes.gtf > tr2g.txt
$ t2g make --version -p - < genes.gtf > tr2g.txt       # (with version number)
```
