/**
 * @class PndMasterSimTask
 * @brief The default sim tasks
 * @details # Master Sim Task Class
 * This class includes all the tasks which need to be used after simulation in the default sim macros. 
 * @remark If you find some obsolete task which needs to be changed, contact the computing coordinator.
 * @author Stefano Spataro <spataro@to.infn.it>, Torino University
 * @version 1.0
 * @date Feb 4, 2016
 **
 **/

#ifndef PNDMASTERSIMTASK_H
#define PNDMASTERSIMTASK_H

#include "PndMasterTask.h"

class TClonesArray;

class PndMasterSimTask : public PndMasterTask
{
 public:

  /** 
   * @brief Default constructor 
   */
  PndMasterSimTask();
  
  /**
   * @brief Destructor
   */
  virtual ~PndMasterSimTask();

  /** 
   * @brief Set the persistency of all the tasks 
   * @param pers Persistency level: 0 no TCA, 1 all TCA
   */
  virtual void SetPersistency(Bool_t pers = kTRUE);
  
 private:
  
  /**
   * @brief enum of the task list
   * @detail This enum avoids to retrieve tasks using the integer value (misleading), but uses an easier scheme. The enum is "k" + the class name, i.e. "kPndEmcHitProducer" for the class PndSimCorrelator. The function PrintTaskList() can be used to check the list of the tasks and the corresponding number.
   */
  enum simTaskList {
    kPndEmcHitProducer
  };
  
  /** @cond CLASSIMP */
  ClassDef(PndMasterSimTask,1);
  /** @endcond */
};

#endif /* PNDMASTERSIMTASK_H */
