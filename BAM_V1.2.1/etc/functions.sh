#! /bin/bash -x

#-------------------------------------------------------------------------------------------------#
#                         Brazilian global Atmospheric Model - BAM_V1.2.1                         #
#-------------------------------------------------------------------------------------------------#
# Descrição:                                                                                      #
#     Script onde são definidas as função para a configuração do modelo BAM                       #
#                                                                                                 #
# Observações:                                                                                    #
#     * Não remova e não mova de lugar as linhas com "#DESCRIPTION"                               #
#     * Ao adicionar novas função, siga o exemplo das funções pré-existentes                      #
#                                                                                                 #
# Revisões:                                                                                       #
#     * 20 Dec 2017 - J. G. de Mattos - versão inicial                                            #
#     * 11 Nov 2019 - C. F. Bastarz - revisão geral                                               #
#                                                                                                 #
# DMD/CPTEC/INPE, 2019                                                                            #
#-------------------------------------------------------------------------------------------------#

assing(){
#DESCRIPTION: Lista as funções (uso interno).
  eval export $1=$2

}

# Variaveis principais
vars_export(){
#DESCRIPTION: Exporta as variáveis de ambiente do BAM.

  FilePaths=$(dirname ${BASH_SOURCE})/paths.conf
  while read line
  do

    assing $line

  done < <(sed -r 's/^\s*(.*\S)*\s*$/\1/;/^$/d;/#/d' ${FilePaths})
}

# Funcao configurar
configurar(){
#DESCRIPTION: Configura o BAM (cria diretórios, links e altera Makefiles).

  # Exportando variaveis
  vars_export

  echo ""
  echo -e "\033[34;1m > A variavel NameInstall possui o valor \033[36;1m${NameInstall}\033[m \033[m"
  echo -e "\033[34;1m > o configurador do BAM ira criar pastas com o nome \033[36;1m${NameInstall}\033[m \033[m"
  echo -e "\033[34;1m > nos discos scratch[in,out] abaixo. \033[m"
  echo ""
  echo -e "\033[34;1m > BAM HOME: \033[36;1m${home_bam}\033[m \033[m"
  echo -e "\033[34;1m > BAM SUBM: \033[36;1m${subt_bam}\033[m \033[m"
  echo -e "\033[34;1m > BAM WORK: \033[36;1m${work_bam}\033[m \033[m" 

  if [ "/"${home_bam}"/" == "/"${HOME}/${NameInstall}"/" ] 
  then 
     echo ""
     echo -e "\033[31;1m > O sistema detectou que voce esta tentando instalar o \033[m"
     echo -e "\033[31;1m > BAM no home, o que nao e recomendado, pois podera haver falta de espaco. \033[m"
  fi

  echo ""

  read -p "Deseja continuar? (S/N) " -n 1 -r

  echo ""

  if [[ $REPLY = [Ss] ]]
  then

    # Mensagem ao usuario
    echo ""
    echo -e "\033[34;1m > Configurando o BAM... \033[m"
  
    # Criando estrutura de diretorios
  
    # Cria datainout em diante 
    if test ! -s ${subt_dataout}; then mkdir -p ${subt_dataout}; fi
    if test ! -s ${work_dataout}; then mkdir -p ${work_dataout}; fi
  
    if test ! -e ${subt_bam}; then mkdir -p ${subt_bam}; fi
    if test ! -s ${work_bam}; then mkdir -p ${work_bam}; fi
  
    # Cria datarun em diante 
    if test ! -s ${subt_run}; then mkdir -p ${subt_run}; fi
    if test ! -s ${work_run}; then mkdir -p ${work_run}; fi
  
    # Cria pastas do BAM
    if test ! -s ${work_pre_bam}; then mkdir -p ${work_pre_bam}; fi
    if test ! -s ${subt_pre_bam_datain}; then mkdir -p ${subt_pre_bam_datain}; fi
    if test ! -e ${work_pre_bam_datain}; then ln -sf ${subt_pre_bam_datain} ${work_pre_bam_datain}; fi
    if test ! -s ${work_pre_bam_dataout}; then mkdir -p ${work_pre_bam_dataout}; fi
    if test ! -e ${subt_pre_bam_dataout}; then ln -sf ${work_pre_bam_dataout} ${subt_pre_bam_dataout}; fi
  
    if test ! -s ${subt_pre_bam_datasst}; then mkdir -p ${subt_pre_bam_datasst}; fi
  
    if test ! -s ${subt_pre_bam_run}; then mkdir -p ${subt_pre_bam_run}; fi
    if test ! -e ${work_pre_bam_run}; then ln -sf ${subt_pre_bam_run} ${work_pre_bam_run}; fi
    if test ! -s ${subt_pre_bam_databcs}; then mkdir -p ${subt_pre_bam_databcs}; fi
    if test ! -e ${work_pre_bam_databcs}; then ln -sf ${subt_pre_bam_databcs} ${work_pre_bam_databcs}; fi
  
    if test ! -s ${work_model_bam}; then mkdir -p ${work_model_bam}; fi
    if test ! -s ${subt_model_bam_datain}; then mkdir -p ${subt_model_bam_datain}; fi
    if test ! -e ${work_model_bam_datain}; then ln -sf ${subt_model_bam_datain} ${work_model_bam_datain}; fi
    if test ! -s ${work_model_bam_dataout}; then mkdir -p ${work_model_bam_dataout}; fi
    if test ! -e ${subt_model_bam_dataout}; then ln -sf ${work_model_bam_dataout} ${subt_model_bam_dataout}; fi
    if test ! -s ${subt_model_bam_run}; then mkdir -p ${subt_model_bam_run}; fi
    if test ! -e ${work_model_bam_run}; then ln -sf ${subt_model_bam_run} ${work_model_bam_run}; fi
  
    if test ! -s ${work_pos_bam}; then mkdir -p ${work_pos_bam}; fi
    if test ! -s ${subt_pos_bam_datain}; then mkdir -p ${subt_pos_bam_datain}; fi
    if test ! -e ${work_pos_bam_datain}; then ln -sf ${subt_pos_bam_datain} ${work_pos_bam_datain}; fi
    if test ! -s ${work_pos_bam_dataout}; then mkdir -p ${work_pos_bam_dataout}; fi
    if test ! -e ${subt_pos_bam_dataout}; then ln -sf ${work_pos_bam_dataout} ${subt_pos_bam_dataout}; fi
    if test ! -s ${subt_pos_bam_run}; then mkdir -p ${subt_pos_bam_run}; fi
    if test ! -e ${work_pos_bam_run}; then ln -sf ${subt_pos_bam_run} ${work_pos_bam_run}; fi
  
    if test ! -s ${subt_grh_bam}; then mkdir -p ${subt_grh_bam}; fi
    if test ! -s ${work_grh_bam_dataout}; then mkdir -p ${work_grh_bam_dataout}; fi
    if test ! -e ${subt_grh_bam_dataout}; then ln -sf ${work_grh_bam_dataout} ${subt_grh_bam_dataout}; fi
    if test ! -s ${subt_grh_bam_run}; then mkdir -p ${subt_grh_bam_run}; fi
    if test ! -e ${work_grh_bam_run}; then ln -sf ${subt_grh_bam_run} ${work_grh_bam_run}; fi
  
    # Altera scripts e Makefiles do BAM
  
    sed -i "/# Caminho onde deverao ficar os arquivos executaveis do Pre/{n;d}" ${home_pre_bam_source}/Makefile
    sed -i "/# Caminho onde deverao ficar os arquivos executaveis do Pre/a\PATH2="${home_pre_bam_run}"" ${home_pre_bam_source}/Makefile
  
    sed -i "/# Caminho onde devera fica o arquivo executavel do Model/{n;d}" ${home_model_bam_source}/Makefile
    sed -i "/# Caminho onde devera fica o arquivo executavel do Model/a\PATH2="${home_model_bam_run}"" ${home_model_bam_source}/Makefile
  
    sed -i "/# Caminho onde deverao ficar os arquivos executaveis do Pos/{n;d}" ${home_pos_bam_source}/Makefile
    sed -i "/# Caminho onde deverao ficar os arquivos executaveis do Pos/a\PATH2="${home_pos_bam_run}"" ${home_pos_bam_source}/Makefile
  
    sed -i "/# Caminho onde devera fica o arquivo executavel do grid history/{n;d}" ${home_grh_bam_source}/Makefile
    sed -i "/# Caminho onde devera fica o arquivo executavel do grid history/a\PATH2="${home_grh_bam_run}"" ${home_grh_bam_source}/Makefile
  
    sed -i "/DirInPut='\//,1d" ${home_run_bam}/PostGridHistory.nml
    sed -i "/!DirInPut=subt_model_bam_dataout/a\DirInPut=\'"${subt_model_bam_dataout}"\/TQ0042L028\'\," ${home_run_bam}/PostGridHistory.nml
  
    sed -i "/DirOutPut='\//,1d" ${home_run_bam}/PostGridHistory.nml
    sed -i "/!DirOutPut=subt_grh_bam_dataout/a\DirOutPut=\'"${subt_grh_bam_dataout}"\/TQ0042L028\'\," ${home_run_bam}/PostGridHistory.nml
  
    sed -i "/DirMain='\//,1d" ${home_run_bam}/PostGridHistory.nml
    sed -i "/!DirMain=subt_bam/a\DirMain=\'"${subt_bam}"\'\," ${home_run_bam}/PostGridHistory.nml
  
    sed -i "/export PATHBASE=\//,1d" ${home_run_bam}/EnvironmentalVariablesMCGA
    sed -i "/# Caminho do BAM no HOME/a\export PATHBASE="${home_bam}"" ${home_run_bam}/EnvironmentalVariablesMCGA
  
    sed -i "/export DK=\//,1d" ${home_run_bam}/EnvironmentalVariablesMCGA
    sed -i "/# Caminho do BAM no scratchin (SUBMIT_HOME)/a\export DK="${subt_bam}"" ${home_run_bam}/EnvironmentalVariablesMCGA
  
    sed -i "/export DK2=\//,1d" ${home_run_bam}/EnvironmentalVariablesMCGA
    sed -i "/# Caminho do BAM no scratchout (SUBMIT_WORK)/a\export DK2="${subt_bam}"" ${home_run_bam}/EnvironmentalVariablesMCGA
  
    echo ""
    echo -e "\033[34;1m > Configuracao BAM completa! \033[m"
  
  elif [[ ${resposta} = N || ${resposta} = "n" ]]
  then
  
    echo -e "\033[34;1m > Saindo do configurador. \033[m"
    exit 0
  
  else
  
    exit 1
  
  fi

}

# Funcao compilar
compilar(){
#DESCRIPTION: Compila o sistema completo (cria executáveis do pré, model e pós).

  # Exportando variaveis principais
   vars_export

  # Verificando se esta logado no eslogin01
  if [ ${HOSTNAME} != "eslogin01" -a ${HOSTNAME} != "eslogin02" ]
  then

    echo ""

    echo "###########################################################################"
    echo "#                                                                         #"
    echo "#                    Você está logado na ${HOSTNAME}                        #"
    echo "#                                                                         #"
    echo "# Antes de proceder com a instalação, você deve logar em um destes hosts: #"
    echo "#                                                                         #"
    echo "# $ ssh eslogin01 -XC                                                     #"
    echo "#                                                                         #"
    echo "#  ou                                                                     #"
    echo "#                                                                         #"
    echo "# $ ssh eslogin02 -XC                                                     #"
    echo "#                                                                         #"
    echo "###########################################################################"
    
    echo ""

    exit 1

  fi
 
  #############################################################################
  # Compilação do BAM
  #############################################################################

  echo ""
  echo "%%%%%%%%%%%%%%%%%%%%"
  echo " Compilando o BAM:  "
  echo "%%%%%%%%%%%%%%%%%%%%"
  echo ""

  #########################################
  # Compilação pre/pos
  #########################################

  echo "+++++++++++++++++++++++++"
  echo "   Compilando o pre:     "
  echo "+++++++++++++++++++++++++" 
 
  cd ${home_pre_bam_source}
  make clean pgi_cray
   
  echo "+++++++++++++++++++++++++"
  echo "   Compilando o pos:     "
  echo "+++++++++++++++++++++++++"   
   
  cd ${home_pos_bam_source}
  make clean pgi_cray
  
  #########################################
  # Compilação model
  #########################################

  echo "+++++++++++++++++++++++++"
  echo "   Compilando o model:   "
  echo "+++++++++++++++++++++++++"   

  cd ${home_model_bam_source}
  make clean pgi_cray

  echo -e "\033[34;1m ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \033[m"
  echo -e "\033[34;1m > Compilacao BAM completa.# Verifique possiveis erros no arquivo de log, caso o tenha criado. \033[m"
  echo -e "\033[34;1m ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \033[m"

  #########################################
  # Copiando executaveis
  #########################################

  #
  # exec BAM
  #
  cp -pfrv ${home_pre_bam_run}/* ${subt_pre_bam_run}
  cp -pfrv ${home_model_bam_run}/* ${subt_model_bam_run}
  cp -pfrv ${home_pos_bam_run}/* ${subt_pos_bam_run}

  echo -e "\033[34;1m > Compilacao BAM completa. Verifique possiveis erros no arquivo de log, caso o tenha criado. \033[m"

  banner

}

# Funcao aloca tescase
testcase(){
#DESCRIPTION: Descompacta os arquivos do testcase.

  vars_export

  #
  # copiando arquivos fixos (Já devem estar pois são copiados durante a configuracão)
  #
  
  cp -pvfr ${public_bam}/model/datain/* ${subt_model_bam_datain}/

  cp -pvfr ${public_bam}/pre/datain/*201[2-3]* ${subt_pre_bam_datain}/
#  cp -pvfr ${public_bam}/pre/datain/*2020* ${subt_pre_bam_datain}/

  cp -pvfr ${public_bam}/pre/dataout/* ${subt_pre_bam_dataout}/
  cp -pvfr ${public_bam}/pre/datasst   ${subt_pre_bam}/
  cp -pvfr ${public_bam}/pre/databcs   ${subt_pre_bam}/
  cp -pvfr ${public_bam}/pre/dataco2   ${subt_pre_bam}/

}

# Final
banner(){
#DESCRIPTION: Imprime algumas informações sobre o script (uso interno).

  echo -e ""
  echo -e "\e[34;1m > Para mais informações sobre esta distribuição do BAM, leia o arquivo: \e[m"
  echo -e "\e[32;1m > ${PWD}/README \e[m"
  echo -e ""

}

ajuda(){
#DESCRIPTION: Mostra a este menu de ajuda.

  echo ""
  echo " Uso.....: ./${0##*/} <opcao>"
  echo ""

  First=1

  grep -i '(){$' ${BASH_SOURCE} | sed 's/(){//g' | while read function; do

    dsc=$(sed -n "/${function}(){$/ {n;p}" ${BASH_SOURCE} | sed "s/#*[DdEeSsCcRrIiPpTtIiOoNn]*://g")

    if [ ${First} -eq 1 ]
    then

        First=0
        printf " Opções..:%2s* \e[1;31m%s\e[m \e[1;37;1m-->\e[m\e[1;34m%s\e[m\n" " " "$function" "${dsc}"

     else

        printf "%12s* \e[1;31m%s\e[m \e[1;37;1m-->\e[m\e[1;34m%s\e[m\n" " " "$function" "${dsc}"

     fi

  done

  echo ""

  First=1

  grep -i '(){$' ${BASH_SOURCE} | sed 's/(){//g' | while read function; do

    dsc=$(sed -n "/${function}(){$/ {n;p}" ${BASH_SOURCE} | sed "s/#*[DdEeSsCcRrIiPpTtIiOoNn]*://g")

    if [ ${First} -eq 1 ]
    then

        First=0
        printf " Exemplos:%2s\e[1;37m${0##*/}\e[m \e[1;31m%s\e[m\n" " " "$function"

     else

        printf "%12s\e[1;37m${0##*/}\e[m \e[1;31m%s\e[m\n" " " "$function"

     fi

  done

}

teste(){
#DESCRIPTION: Função para testar uso de funções.
  echo ""
  echo " Olá mundo!"

}
