#! /bin/bash

#-------------------------------------------------------------------------------------------------#
#                         Brazilian global Atmospheric Model - BAM_V1.2.1                         #
#-------------------------------------------------------------------------------------------------#
# Descrição:                                                                                      #
#     Script para a configuração e a compilação do modelo BAM                                     #
#                                                                                                 #
# Uso:                                                                                            #
#     ./config_smg.sh <opção>                                                                     #
#                                                                                                 #
# Opções:                                                                                         #
#     * assing..........: lista as funções (uso interno)                                          #
#     * vars_export.....: exporta as variáveis de ambiente do BAM                                 #
#     * copy_fixed_files: copia os arquivos fixos necessários para qualquer rodada                #
#     * configurar......: configura o BAM (cria diretórios, links e altera Makefiles)             #
#     * compilar........: compila o sistema completo (cria executáveis do pré, model e pós)       #
#     * testcase........: descompacta os arquivos do testcase                                     #
#     * banner..........: imprime algumas informações sobre o script (uso interno)                #
#     * ajuda...........: mostra a este menu de ajuda                                             #
#     * teste...........: função para testar uso de funções                                       #
#                                                                                                 #
# Exemplos:                                                                                       #
#     * Configura a estrutura de diretórios do BAM                                                #
#       ./config_bam.sh configurar                                                                #
#     * Compila o pré-processamento, modelo e pós-processamento (PGI, utilizado apenas na Tupã)   #
#       ./config_bam.sh compilar                                                                  #
#     * Copia e aloca os dados do testcase (utilizado apenas na Tupã)                             #
#       ./config_bam.sh testcase                                                                  #
#                                                                                                 #
# Observações:                                                                                    #
#     * Os caminhos do BAM estao contidos no arquivo etc/paths.sh                                 #
#     * As funcoes que podem ser executadas estao no arquivo etc/functions.sh                     #
#                                                                                                 #
# Revisões:                                                                                       #
#     * 20 Dec 2017 - J. G. de Mattos - separado em diferentes arquivos                           #
#                                       * etc/paths.shi: caminhos                                 #
#                                       * etc/functions.sh: funcoes                               #
#     * 03 Dec 2017 - J. G. de Mattos - Adaptando para a instalação                               #
#                                       standalone do bam                                         #
#                                                                                                 #
#     * 11 Nov 2019 - C. F. Bastarz   - ajustes para a instalação e configuração                  #
#                                       do BAM.                                                   #
#                                                                                                 #
# DMD/CPTEC/INPE, 2019                                                                            #
#-------------------------------------------------------------------------------------------------#

RootDir=$(dirname ${BASH_SOURCE})
. ${RootDir}/etc/functions.sh

#
# Verifica os argumentos passados junto com o script
# 

echo -e ""
echo -e "\e[36;1m >>> ${BASH_SOURCE##*/} executado a partir de\e[m \e[32;1m${0##*/}\e[m"

if [ $# = 0 ];then
  echo -e ""
  echo -e "\e[31;1m > Nao foi passado nenhum argumento! \e[m"
  ajuda
  banner
  exit -1
fi

echo -en "\e[34;1m Opcao escolhida: \e[m \e[37;1m ${1} \e[m"

f=0
for function in $(grep -i '(){$' ${RootDir}/etc/functions.sh | sed 's/(){//g');do
   if [ ${1} == ${function} ];then
     f=1
     echo -e "\e[37;1m[\e[m\e[32;1m OK \e[m\e[37;1m]\e[m"
     vars_export
     ${1}
     banner
   fi
done

if [ ${f} -eq 0 ];then
  
  echo -e "\e[37;1m[\e[m\e[31;1m FAIL \e[m\e[37;1m]\e[m"
  echo -e ""
  echo -e "\e[37;1m Opcao desconhecida, <ajuda>: \e[m"
  vars_export
  ajuda
  banner
fi
