#include "processingTool_getInteractionGroups.hpp"

void getInteractionGroups(Rockable *box, std::vector<size_t>& nbInt) {
  nbInt.clear();
  std::sort(box->activeInteractions.begin(), box->activeInteractions.end(), std::less<Interaction*>());
  size_t igrp = box->activeInteractions[0]->i;
  size_t jgrp = box->activeInteractions[0]->j;
  size_t nb = 1;
  for (size_t k = 1; k < box->activeInteractions.size(); k++) {
    if (box->activeInteractions[k]->i == igrp && box->activeInteractions[k]->j == jgrp) {
      nb++;
    } else {
      nbInt.push_back(nb);
      igrp = box->activeInteractions[k]->i;
      jgrp = box->activeInteractions[k]->j;
      nb = 1;
    }
  }
  nbInt.push_back(nb);
}