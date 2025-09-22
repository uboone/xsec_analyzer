SAMDEF=prodnuwro_nu_overlay_run2_pandora_reco2_run2b_reco2
( samweb list-definition-files -e uboone ${SAMDEF} | while read FILENAME; do PREFIX='enstore:'; SUFFIX='(.*)'; DIR_TO_FILE=$(samweb -e uboone locate-file ${FILENAME} | sed -e "s@${PREFIX}@@" -e "s@${SUFFIX}@@"); echo ${DIR_TO_FILE}/${FILENAME}; done ) 2>&1 | tee /pnfs/uboone/persistent/users/birwin/myLists/run2b_nuwro_fake_data.list
SAMDEF=prodnuwro_nu_overlay_run3_pandora_reco2_run3a_reco2
( samweb list-definition-files -e uboone ${SAMDEF} | while read FILENAME; do PREFIX='enstore:'; SUFFIX='(.*)'; DIR_TO_FILE=$(samweb -e uboone locate-file ${FILENAME} | sed -e "s@${PREFIX}@@" -e "s@${SUFFIX}@@"); echo ${DIR_TO_FILE}/${FILENAME}; done ) 2>&1 | tee /pnfs/uboone/persistent/users/birwin/myLists/run3a_nuwro_fake_data.list
SAMDEF=prodnuwro_nu_overlay_run3_pandora_reco2_run3b_reco2
( samweb list-definition-files -e uboone ${SAMDEF} | while read FILENAME; do PREFIX='enstore:'; SUFFIX='(.*)'; DIR_TO_FILE=$(samweb -e uboone locate-file ${FILENAME} | sed -e "s@${PREFIX}@@" -e "s@${SUFFIX}@@"); echo ${DIR_TO_FILE}/${FILENAME}; done ) 2>&1 | tee /pnfs/uboone/persistent/users/birwin/myLists/run3b_nuwro_fake_data.list
SAMDEF=nuwro_overlay_bnb_reco2_pandora_run4_reco2
( samweb list-definition-files -e uboone ${SAMDEF} | while read FILENAME; do PREFIX='enstore:'; SUFFIX='(.*)'; DIR_TO_FILE=$(samweb -e uboone locate-file ${FILENAME} | sed -e "s@${PREFIX}@@" -e "s@${SUFFIX}@@"); echo ${DIR_TO_FILE}/${FILENAME}; done ) 2>&1 | tee /pnfs/uboone/persistent/users/birwin/myLists/run4_nuwro_fake_data.list
