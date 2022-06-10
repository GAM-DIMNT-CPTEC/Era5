# Compilação do BAM_V1.2.1 na EGEON

## Ajuste dos Makefile

Para a compilação do pre, model e pos, foi aproveitado o arquivo `Makefile.linux_gnu`, com as seguintes modificações:

```
FC = mpif90
```

Este arquivo com as modificações, foi renomeado para `Makefile.linux_gnu_egeon` e para o model, foi removida a opção `-static` para a compilação, além de adicionar a opção `-ffree-line-length-none`.

## Compilação do pre

Para compilar o pre, é necessário compilar a `libbacio`. Para isso, foi necessário modificar o script `makebacio_cray_gnu.sh`, incluindo as instruções a seguir:

```
gnu)
  export FCMP=${1:-mpif90}
  export CCMP=${2:-cc}
  flagOpt="-O0"
  flag64bit=""
```

Este script foi renomeado para `makebacio_cray_gnu-carlos.sh` e foi utilizado para a compilação. É importante notar que este script não é utilizado pelo arquivo `Makefile.linux_gnu_egeon` e que a compilação do pre, deve ser feita sem a opção `clean`, como em:

```
make linux_gnu_egeon
```

## Ambiente de compilação para o pre e pos

```
Currently Loaded Modules:
  1) gnu9/9.4.0   3) openmpi4/4.1.1   5) netcdf-fortran/4.5.3   7) hwloc/2.5.0
  2) ucx/1.11.2   4) netcdf/4.7.4     6) phdf5/1.10.8           8) libfabric/1.13.0
```

## Compilação do pos

Basta compilar normalmente com o comando:

```
make clean linux_gnu_egeon
```

## Compilação do model

O model foi compilado com uma versão mais antiga do GNU da EGEON e foi utilizado o container `ubuntu_remix_latest-gcc_4.8.5.sif`. Isso foi feito porque com a versão do GNU da EGEON, a compilação apresenta problemas e para apenas testar o método de perturbação do oensMB09, foi suficiente aproximar o ambiente de compilação utilizado na Tupã (o qual é fornecido pelo container).

Para compilar:

```
module load singularity
singularity shell -e --bind /mnt/beegfs/carlos.bastarz:/mnt/beegfs/carlos.bastarz ubuntu_remix_latest-gcc_4.8.5.sif
cd BAM_V1.2.1/model/source
make clean linux_gnu_egeon
```
