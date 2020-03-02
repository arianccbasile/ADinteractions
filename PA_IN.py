
# coding: utf-8

# In[1]:


from os.path import expanduser, join
import os
import cheroot
import pkg_resources
import json
import cherrypy
import cobra
import os
import pandas as pd
from multiprocessing import Process, Pool
from itertools import combinations

# In[2]:
def get_all_pairs(source_models):
    """ Get all of the unique pairs from a list of models.
    Parameters
    ----------
    source_models : list of str
        List of path names to model files
    Returns
    -------
    list of tuple
        List of tuples where each tuple has two elements, path to first model file in pair and
        path to second model file in pair
    """

    return [pair for pair in combinations(source_models, 2)]


'''
    The set of functions in this widget allow the automatic creation of two-species community metabolic models, either from user defined list of pairs of species or from a list created automatically from the models available in the folder with individuals species metabolic models.
'''


def totalEXRxns(modelA,modelB): 
    '''
    This function creates a list object with the id (Reaction.id method from cobrapy) for all the unique exchange reactions found in both models listed (modelA, modelB). These exchange reactions are then differentiated from the reactions from models A and B by an additional tag in the end of the reaction id ([u])
    :param modelA: cobrapy Model object
    :param modelB: cobrapy Model object
    :return EX_finalRxns: list with the ids of all the unique exchange reactions from modelA and modelB.
    '''

    #cherrypy.log('Started function totalEXRxns. Working on %s as modelA and %s as modelB' %(modelA,modelB))


    # List all the exchange reactions in modelA
    EX_rxnsA = set()
    for i in range(len(modelA.reactions)):
        rxnsA = str(modelA.reactions[i])
        if 'EX_' in rxnsA:
            EX_rxnsA.add(rxnsA)
    

    #cherrypy.log('Finished finding all exchange reactions in modelA. There are %d of them' %(len(EX_rxnsA)))


    # List all the exchange reactions in modelB
    EX_rxnsB = set()
    
    for j in range(len(modelB.reactions)):
        rxnsB = str(modelB.reactions[j])
        if 'EX_' in rxnsB:
            EX_rxnsB.add(rxnsB)


    #cherrypy.log('Finished finding all exchange reactions in modelB. There are %d of them' %(len(EX_rxnsB)))


    # List all the different exchange reactions that are present in models A and B. They will have some that overlap.
    EX_total =  list(EX_rxnsA | EX_rxnsB)


    #cherrypy.log('Finished creating a list with all the exchange reactions existing in the two models. There are %d of them' %(len(EX_total)))


    # Create a list with all the Exchange reactions and assign them a new identifier. These will be the exchange reactions that will make up the external compartment that is common to both bacterial species.
    EX_finalRxns = []
    for each in range(len(EX_total)):
        rxn = EX_total[each] + '[u]'
        EX_finalRxns.append(rxn)


    #cherrypy.log('Finished adding the tag [u] to the end of the exchange reactions')

    return EX_finalRxns


def createEXmodel(EXreactions): 
    '''
    This function takes the list of exchange reactions created using the function totalEXRxns and creates a Model object using cobrapy composed of those reactions with the upper bound flux values of 1000, lower bound flux values of -1000, and objective coefficient of 0, and one metabolite as being uptaken by the reaction (stoichiometric coefficient of -1). This is a model composed solely of exchange reactions and it's the model for the extra compartment created for the full community model
    :param EXreactions: list of reactions that are the output of the function totalEXRxns (above)
    :return exchange_model: cobrapy Model object for the compartment that will serve as an extra compartment in the full community model.
    '''


    #cherrypy.log("Started the function that creates the exchange reactions for the community model")

    exchange_model = cobra.Model('Model with the exchange reactions only')

    #cherrypy.log("Created the base exchange model object")
    
    for i in EXreactions: 
        new_i = str(i)
        new_i = new_i[3:]
        new_met = cobra.Metabolite(new_i)
        
        rxn = cobra.Reaction(i)
        rxn.lower_bound = -1000.000
        rxn.upper_bound = 1000.000
        rxn.objective_coefficient = 0.000
        rxn.add_metabolites({new_met:-1.0}) 
        
        exchange_model.add_reaction(rxn)


    #cherrypy.log('Finished adding all exchange reactions in exchange model object. There are %d of them' %(len(exchange_model.reactions)))

    return exchange_model 


def createReverseEXmodel(EXreactions): 
    '''
    This function takes the list of exchange reactions created using the function totalEXRxns and creates a Model object using cobrapy composed of those reactions with the upper bound flux values of 1000, lower bound flux values of -1000, and objective coefficient of 0, and one metabolite as being produced by the reaction (stoichiometric coefficient of 1). This is a model composed solely of exchange reactions. The metabolite information for these reactions will be used to update the metabolites of the exchange reactions for models A and B.
    :param EXreactions: list of reactions that are the output of the function totalEXRxns (above)
    :return exchange_modelRev: cobrapy Model object containing only exchange reactions with the production of their respective metabolites
    '''

    #cherrypy.log("Started the function that creates the reverse exchange reactions for the community model")

    exchange_modelRev = cobra.Model('Model with the exchange reactions only with reversed stoi coefficient')

    #cherrypy.log("Created the base reverse exchange model object")

    for i in EXreactions: 
        new_i = str(i)
        new_i = new_i[3:]
        new_met = cobra.Metabolite(new_i)
        
        rxn = cobra.Reaction(i)
        rxn.lower_bound = -1000.000
        rxn.upper_bound = 1000.000
        rxn.objective_coefficient = 0.000
        rxn.add_metabolites({new_met:1.0})
        
        exchange_modelRev.add_reaction(rxn)

    #cherrypy.log('Finished adding all exchange reactions in reverse exchange model object. There are %d of them' %(len(exchange_modelRev.reactions)))

    return exchange_modelRev



def addEXMets2SpeciesEX(reverseEXmodel,speciesModel):
    '''
    This function takes the model with exchange reactions where the metabolite is produced (output from function createReverseEXmodel) and a species model, and adds the metabolite from the reverse model to the exhange reactions of the species model. For instance:
    Reaction :  modelB_EX_cpd11588_e0 got the cpd11588_e0[u] added.
                'model_B_cpd11588_e0 <=> cpd11588_e0[u]'
    This way, when a compound is exported to the extracellular environment, it is automatically transformed into a form that is common to all members in the community.
    :param reverseEXmodel: cobrapy Model object containing only exchange reactions with the production of their respective metabolites
    :param speciesModel: Model object of a particular species.
    :return speciesModel: Model object of a particular species with updated exchange reactions are updated.
    '''

    #cherrypy.log('Started function to add metabolites to the exchange reactions of the reverse exchange model') #not right

    for j in range(len(reverseEXmodel.reactions)):
        exRxn = str(reverseEXmodel.reactions[j])
        
        for i in range(len(speciesModel.reactions)):
            rxn = str(speciesModel.reactions[i])
            if rxn in exRxn:
                new_met = reverseEXmodel.reactions[j].metabolites 
                speciesModel.reactions[i].add_metabolites(new_met)
                speciesModel.reactions[i].lower_bound = -1000.000
                speciesModel.reactions[i].upper_bound = 1000.000

    #cherrypy.log('Finished adding metabolites to the exchange reactions of the reverse exchange model')
    return speciesModel       
               

def replaceRxns(model,modelID):    
    '''
    This function adds the tag specified in the parameter modelID to the beginning of the reaction IDs for a particular Model object. We are doing this so that we know which reactions come from one species or the other. This is the same as assigning each species to a different compartment. This is important because the two species have common reactions and metabolites, but are not sharing these metabolites in their biology, since the cells are closed compartments. They only share the metabolites that are transported in and out of the cell, hence the creation of an extra external compartment.
    :param model: Model object containing the metabolic model of a particular species
    :param modelID: Tag to add to the beginning of the reaction IDs of the model.
    :return model: same model but with updated reactions IDs
    '''


    #cherrypy.log('Started function to replace the reaction IDs in the species models')

    
    for i in range(len(model.reactions)):
        old_rxns = str(model.reactions[i])
        new_rxns = 'model' + modelID + '_' + old_rxns
        model.reactions[i].id = new_rxns

    #cherrypy.log('Finished changing the reaction IDs in the species models')

def replaceMets(model,modelID):
    '''
    This function adds the tag specified in the parameter modelID to the beginning of the metabolite IDs for a particular Model object. We are doing this so that we know which metabolites come from one species or the other. This is the same as assigning each species to a different compartment. This is important because the two species have common reactions and metabolites, but are not sharing these metabolites in their biology, since the cells are closed compartments. They only share the metabolites that are transported in and out of the cell, hence the creation of an extra external compartment.
    :param model: Model object containing the metabolic model of a particular species
    :param modelID: Tag to add to the beginning of the metabolite IDs of the model.
    :return model: same model but with updated metabolite IDs
    '''

    #cherrypy.log('Started function to replace the metabolite IDs in the species models')

    
    for i in range(len(model.metabolites)):
        old_mets = str(model.metabolites[i])
        new_mets = 'model_' + modelID + '_' + old_mets
        
        model.metabolites[i].id = new_mets

    #cherrypy.log('Finished changing the metabolite IDs in the species models')



def createCommunityModel(modelFileA, modelFileB, comFolder):
    '''
    This function takes advantage of the outputs of all the functions defined previously to actually piece together the individual species models and the extra compartment model.
    :param modelFileA: path to the metabolic model of species A in SBML format
    :param modelFileB: path to the metabolic model of species B in SBML format
    :param comFolder:  path to the folder where the metabolic models of the two-species communities will be stored.
    :return two-species community model: in SBML format exported to the folder designated by the user (comFolder) to store these models
    '''
    
    #cherrypy.log('Started function to create community models. ModelFileA is %s, ModelFileB is %s, and the folder where we are going to put the files in is %s' %(modelFileA, modelFileB, comFolder))

    
    #import the model file into the Model object model1 using cobrapy
    try:
        if modelFileA.endswith('.mat'):
            #cherrypy.log('The extension is .mat .This is modelfileA %s' %modelFileA)
            model1 = cobra.io.load_matlab_model(modelFileA)
        elif modelFileA.endswith('.xml') or modelFileA.endswith('.sbml'):
            #cherrypy.log('The extension is .xml or .sbml .This is modelfileA %s' %modelFileA)
            model1 = cobra.io.read_sbml_model(modelFileA)
        elif modelFileA.endswith('.json'):
            #cherrypy.log('The extension is .json . This is modelfileA %s' %modelFileA)
            model1 = cobra.io.load_json_model(modelFileA)
        else:
            #cherrypy.log('We were not able to find a model. This is modelfileA %s' %modelFileA)
            print "not able to find model %s" %modelFileA
    except Exception as e:
        print e

    #cherrypy.log('%s loaded successfully' %model1.id)

    #import the model file into the Model object model2 using cobrapy
    if modelFileB.endswith('.mat'):
        #cherrypy.log('The extension is .mat .This is modelfileB %s' %modelFileB)
        model2 = cobra.io.load_matlab_model(modelFileB)
    elif modelFileB.endswith('.xml') or modelFileB.endswith('.sbml'):
        #cherrypy.log('The extension is .xml or .sbml .This is modelfileB %s' %modelFileB)
        model2 = cobra.io.read_sbml_model(modelFileB)
    elif modelFileB.endswith('.json'):
        #cherrypy.log('The extension is .json . This is modelfileB %s' %modelFileB)
        model2 = cobra.io.load_json_model(modelFileB)
    else:
        #cherrypy.log('We were not able to find a model. This is modelfileB %s' %modelFileB)
        print "not able to find model %s" %modelFileB
        
    #cherrypy.log('%s loaded successfully' %model2.id)


    # Create a communityID to identify the output files belonging to each 2-species community created
    communityID = model1.id+ 'X' + model2.id

    #cherrypy.log('The communityID (reflected in the filename is %s .'%communityID)
    

    # Get all the reactions identified as exchange reactions in both models you're mixing and create a list with of exchange reactions. Then use this list to create what is called an exchange reaction model. This exchange reaction model will be the equivalent of an outside world model, or the lumen for instance, as in the models in Heinken and Thiele AEM 2015 paper. Later manipulation of this particular model will allow the user to choose the diet under which the communities are growing.
    #cherrypy.log('Lets actually run the function that creates te model with the exchange reactions.')
    exModel = createEXmodel(totalEXRxns(model1, model2))

    

    # Create a model that has the fluxes of the exchange reactions reversed. This is because these reactions will will added specifically to each species, in the model. So we are extending the original species models to have more reactions so that each species can exchange metabolites with exchange reactions model. The exchange reactions models then becomes a comparment shared by all the other species in the community model. This is what will allow us to determine how the species interact when they have to share resources.
    #cherrypy.log('Lets actually run the function that creates te model with the reverse exchange reactions.')
    revEXmodel = createReverseEXmodel(totalEXRxns(model1, model2))

    # Add a tag to the metabolite IDs of modelA and modelB.
    #cherrypy.log('Lets actually run the function that replaces the metabolite IDs for modelA (%s) and modelB (%s).' %(model1,model2))
    replaceMets(model1,'A')
    replaceMets(model2,'B')

    # Add the metabolites of the external model to the exchange reactions of each species.
    #cherrypy.log('Lets add the metabolites to the exchange reactions of the species models')
    new_m1 = addEXMets2SpeciesEX(revEXmodel,model1) 
    new_m2 = addEXMets2SpeciesEX(revEXmodel,model2) 


    # Add a tag to the reaction IDs of modelA and modelB.
    #cherrypy.log('Lets replace the reaction IDs on the species models')
    replaceRxns(new_m1,'A')
    replaceRxns(new_m2,'B')

    #cherrypy.log('Finished replacing the reaction IDs on the species models')
    
    

    # Actually create the community model. All previous steps were changing the models that will be put together in the community so that the reactions and metabolites for each organism can still be distinguished and there is proper compartmentalization of reactions and metabolites. Because you can't really create a model from the sum of other 2, I just created a new model (mix) that is exactly model1. The alternative would have been to create an empty model, then add the reactions and metabolites of model1, model2, and exModel.
    #cherrypy.log('Now, after we prepared every piece, we are going to put it all together in a community model')
    mix = new_m1
    mix.id = communityID
    mix.add_reactions(new_m2.reactions)
    mix.add_metabolites(new_m2.metabolites)
    mix.add_reactions(exModel.reactions)
    mix.add_metabolites(exModel.metabolites)

    #cherrypy.log('A community model with the id %s was created. It has %d reactions and %d metabolites'%(mix.id, len(mix.reactions),len(mix.metabolites)))

    # Export the newly created community model to its folder. The models should then be ready to be further analyzed on Widget 5
    cobra.io.write_sbml_model(mix, "%s/community%s.sbml" %(comFolder,communityID))

    cherrypy.log('The model has been exported to the %s folder'%comFolder)


    

def allPairComModels(listOfPairs,modelFolder,comFolder):
    '''
    This function goes through a list with the models that should be paired together to form a community and creates the corresponding two-species community metabolic model using the function createCommunityModel
    :param listOfPairs: file with pairs of species that will make up each 
    two-species community metabolic model.
    :param modelFolder: path to the folder containing the metabolic models of individual species in a SBML format
    :param comFolder: path to the folder that will store the two-species community metabolic models.
    :return set of two-species community metabolic models
    '''
    import os

    #cherrypy.log('We are going to run the createCommunityModel script for all pairs of models in the %s file' %listOfPairs)

    # Check if the directory where two-species community models will be stored exists. If not, create it.
    if not os.path.exists(comFolder):
        os.makedirs(comFolder)
    
    pairsListFile = open(listOfPairs,'r')
    pairsList = []


    for i in pairsListFile:
        i = i.rstrip()
        i = i.replace("'","")
        i = i.split()
        pairsList.append(i)
    
    #cherrypy.log('We created a list with the list of model pairs that will be put together. This list has %d pairs.' %(len(pairsList)))

    # Go down the list containing pairs of species and create a two-species community model.
    for i in range(len(pairsList)):
        modelA = modelFolder + '%s' %pairsList[i][0]
        modelA = str(modelA)
        #cherrypy.log('The pair number %d will use file %s as modelA.' %(i,modelA))

        modelB = modelFolder + '%s' %pairsList[i][1]
        modelB =str(modelB)
        #cherrypy.log('The pair number %d will use file %s as modelB.' %(i,modelB))
        try:
            createCommunityModel(modelA,modelB,comFolder)
        except Exception as e:
            print e
    
    #cherrypy.log('We finished creating the models for all pairs in your list.')

    pairsListFile.close()


# In[3]:



def getListOfModels(comFolder):
    '''
    This function creates a list with all the community models that will be used in the analysis. It creates this list by listing the metabolic models in SBML format present in the user specified folder that contains the community models.
    :param comFolder: path to the folder that contains the community metabolic models.
    :return listOfModels: list object containing the filenames for all the models that will be analysed in the function calculateGR
    '''

    import os
    #cherrypy.log('We will first get the full list of community models from the %s folder' %comFolder)
    
    path = comFolder+"/"
    listOfFiles = os.listdir(path)

    listOfModels = []
    
    for file in listOfFiles:
        if (file.endswith('.sbml') or file.endswith('.xml')):
            pathToFile = path + file
            listOfModels.append(pathToFile)

    #cherrypy.log('There are %s community models what will be analyzed.'%listOfModels)

    return listOfModels


def calculateGR(diet, comFolder, OutFile="OutputGR.txt"):
    '''
    In this function we use cobrapy to calculate the growth rates of the two species that make up the two species community metabolic models under particular metabolite availability conditions. We start by loading the community model (full model) into 3 distinct Model objects with cobrapy. We then change the fluxes of the exchange reactions of the external model so they have lower bounds corresponding to whichever 'Diet' condition the user specifies. We then run a flux balance analysis on the full model, optimizing the biomass reactions of the two species that make up the community at the same time. It then remove all reactions whose IDs start with modelA from the model in modelMinusA, thus leaving only the reactions from modelB and from the external compartment. It then runs a FBA on it, maximizing the biomass reaction for the model with the tag modelB. It then does the samething but for reactions tagged with modelB on modelMinusB. The optimal flux values for the biomass reactions of each species resulting from optimization in the full model and in each model containing only one species, which correspond to predicted growth rates, are then exported to a table in the tab-delimited text formal to a folder chosen by the user.
    :param diet: the metabolite availability conditions. Default on MMinte is complete, but the user can choose another value ('Variant1 through 10')
    :param comFolder: path to the folder containing all the two-species community metabolic models.
    :return outputGRs: table with growth rate information for each of the species belonging to a two-species community metabolic model in the presence and absence of another species.
    '''
    growth_rate_cutoff = 1e-6
    #cherrypy.log('We will now calculate the growth rates of the two species in a community model in the presence and absence of the other species')
    growthRatesFile = open(OutFile,'a+')
    print>>growthRatesFile, 'ModelName', '\t', 'ObjFuntionSpeciesA', '\t', 'ObjFunctionSpeceisB', '\t', 'GRSpeciesAFull','\t', 'GRSpeciesBFull','\t','GRASolo','\t','GRBSolo'
    

    # Create a list of all the models that will be analysed
    allModels = getListOfModels(comFolder)


    for item in range(len(allModels)):
        
        '''
        @summary: load the model, all versions that will be manipulated in the analysis.
        '''
        try:
        # Import the models with cobrapy
            #cherrypy.log(allModels[item])
            modelFull = cobra.io.read_sbml_model(allModels[item])
            modelMinusA = cobra.io.read_sbml_model(allModels[item])
            modelMinusB = cobra.io.read_sbml_model(allModels[item])

            #cherrypy.log('We successfully loaded the file %s into three different Model objects with ids, %s, %s, and %s. They should all have the same id.'%(allModels[item],modelFull.id,modelMinusA.id,modelMinusB.id))
        
        # Determine what the objective function is. It should be composed of two reactions, the biomass reactions for each of the species that compose the model. Store the biomass reactions in a new variable to be used later.
            ObjKeys = modelFull.objective.keys()
            idObjKeys = ObjKeys[0].id, ObjKeys[1].id
        
        # Open the metabolite conditions file ('Diet')
            dietValues = open(diet,'r')
    
        # Default 'Diet' for MMinte is Complete. After choosing the 'Diet' change the lower bounds for the exchange reactions of the external compartment.

#            cherrypy.log('The diet chosen for this particular run of the functions was %s' %diet)
            #cherrypy.log(str(modelFull))
            for line in dietValues:
                try:
                    new_line = line.rstrip('\n').split('\t')
                    modelFull.reactions.get_by_id(new_line[0]).lower_bound = -float(new_line[1])
                    modelMinusA.reactions.get_by_id(new_line[0]).lower_bound = -float(new_line[1])
                    modelMinusB.reactions.get_by_id(new_line[0]).lower_bound = -float(new_line[1])
                    #cherrypy.log(new_line[0]+"found")
                except:
                    #cherrypy.log(new_line[0])
                    continue


            #cherrypy.log('We finished changing the lower bounds for the fluxes of the exchange reactions in the models to better fit the availability of metabolites for the microbial communities we are simulating the growth of. ')
            dietValues.close()

        # Run FBA on Full model
            modelFull.optimize()
        


        #cherrypy.log('The objective function of the full model, %s, contains the objective functions of the two models that make it up. They are %s.'%(allModels[item], idObjKeys))

        #cherrypy.log('We finished running the FBA on the full model. The status of the solution is %s and the value of the solution found is %f'%(modelFull.solution.status,modelFull.solution.f))

        # Find which reactions are tagged as being from modelA, store them in a list, then create the modelMinusA, that is, remove all reactions that are part of one of the species in the model.
        # Run FBA on that reduced model.

            #cherrypy.log('We are now going to remove all reactions that are tagged as being of species A from the model.')

            listSilentItemsA = []
        
            for item in modelMinusA.reactions:
                item = str(item)
                if item.startswith('modelA_'):
                    listSilentItemsA.append(item)

            #cherrypy.log('We found %d reactions tagged as being of species A'%len(listSilentItemsA))


            for j in listSilentItemsA:
                rxn = j.strip()
                deadRxnA = modelMinusA.reactions.get_by_id(rxn)
                deadRxnA.remove_from_model()

            #cherrypy.log('We finished removing all reactions from the model. The model now has %s reactions and the reaction on the objective function in %s'%(len(modelMinusA.reactions),modelMinusA.objective))
        
            modelMinusA.optimize()

            #cherrypy.log('We finished running the FBA on the model without reactions of species A. The status of the solution is %s and the value of the solution found is %f'%(modelMinusA.solution.status,modelMinusA.solution.f))


        # Find which reactions are tagged as being from modelB, store them in a list, then create the modelMinusB, that is, remove all reactions that are part of one of the species in the model.
        # Run FBA on that reduced model.

            #cherrypy.log('We are now going to remove all reactions that are tagged as being of species B from the model.')

            listSilentItemsB = []
        
            for item in modelMinusB.reactions:
                item = str(item)
                if item.startswith('modelB_'):
                    listSilentItemsB.append(item)

            #cherrypy.log('We found %d reactions tagged as being of species B'%len(listSilentItemsB))
            for j in listSilentItemsB:
                rxn = j.strip()
                deadRxnB = modelMinusB.reactions.get_by_id(rxn)
                deadRxnB.remove_from_model()

            #cherrypy.log('We finished removing all reactions from the model. The model now has %s reactions and the reaction on the objective function in %s'%(len(modelMinusB.reactions),modelMinusB.objective))
        
            modelMinusB.optimize()

            #cherrypy.log('We finished running the FBA on the model without reactions of species A. The status of the solution is %s and the value of the solution found is %f'%(modelMinusB.solution.status,modelMinusB.solution.f))


        # Get the x_dict values (fluxes) for the reactions listed under idObjKeys for all three models.
        # Output them to a file with a table that has the information about the model, the species ID in the model, and the growth rates of each species in the full model and in isolation.


            ObjA = []
            ObjB = []
        
            if idObjKeys[0].startswith('modelA'):
                ObjA = idObjKeys[0]
            else:
                ObjB = idObjKeys[0]
            
        
            if idObjKeys[1].startswith('modelB'):
                ObjB = idObjKeys[1]
            else:
                ObjA = idObjKeys[1]

            #cherrypy.log('We matched the reactions in the objective function to the model they came from. %s was originally from modelA and %s was originally from mobelB' %(ObjA,ObjB))
        
            grAfull = modelFull.solution.x_dict[ObjA]
            grBfull = modelFull.solution.x_dict[ObjB]

            #cherrypy.log("We are going to create a table with the information about the growth of A and B alone and in the presence of B and A respectively.")
        
            if ObjA in modelMinusA.solution.x_dict:
                grAMinusA = modelMinusA.solution.x_dict[ObjA]
            else:
                grAMinusA = 'Solo'
        
        
            if ObjB in modelMinusA.solution.x_dict:
                grBMinusA = modelMinusA.solution.x_dict[ObjB]
            else:
                grBMinusA = 'Solo'
        
        
            if ObjA in modelMinusB.solution.x_dict:
                grAMinusB = modelMinusB.solution.x_dict[ObjA]
            else:
                grAMinusB = 'Solo'
        
        
            if ObjB in modelMinusB.solution.x_dict:
                grBMinusB = modelMinusB.solution.x_dict[ObjB]
            else:
                grBMinusB = 'Solo'

            if grAMinusA != 'Solo' or grBMinusB != 'Solo':
                cherrypy.log('There is a problem with the attribution of growth rate values to their respective species in model %s .'%modelFull.id)


            modelID = modelFull.id
            organisms = modelID.split('X')
        #grAfull,'\t', grBfull,'\t',grAMinusB,'\t',grBMinusA
            # Round very small growth rates to zero.
            if grAfull < growth_rate_cutoff:
                grAfull = 0.
            if grBfull < growth_rate_cutoff:
                grBfull = 0.
            if grAMinusB < growth_rate_cutoff:
                grAMinusB = 0.
            if grBMinusA < growth_rate_cutoff:
                grBMinusA = 0.

            print>> growthRatesFile, modelID, '\t', organisms[0], '\t', organisms[1], '\t', grAfull,'\t', grBfull,'\t',grAMinusB,'\t',grBMinusA
            cherrypy.log("next")
        except: 
            #cherrypy.log("model had problems")
            continue
    #cherrypy.log('We finished calculating the growth rates of the species in isolation and when in the presence of another species and dumped the information to the file: %s' %growthRatesFile)
    growthRatesFile.close()




def calculate_growth_rates_multiproc(diet,comFolder,n_processes=32):
    """ Calculate growth rates for the pairs in a community.
    The medium is a dictionary with an exchange reaction ID as the key and the
    absolute value of bound in direction of metabolite creation as the value
    (i.e. lower bound for `met <--` or upper bound for `met -->`). For example,
    {'EX_h2o': 0.0, 'EX_h2s': 1.0, 'EX_pi': 10.0, ...}
    Parameters
    ----------
    pair_models : list of str
        List of path names to two species community model files
    medium : dict
        Dictionary with exchange reaction ID as key and bound as value
    n_processes: int, optional
        Number of processes in job pool
    Returns
    -------
    pandas.DataFrame
        Results of growth rate calculations
    """
    pool = Pool(n_processes)
    result_list = [pool.apply(calculateGR,args=(diet, comFolder+"/"+i))
    for i in os.listdir(comFolder)]
    pool.close()

# In[5]:



def evaluateInteractions(inGRs, outInter):
    '''
    This function goes over the file with the growth rates of the species that make up a two-species community model and determines the kind of interaction occurring in between the two species. The types interactions are determined according to the paper by Heinken and Thiele 2015 AEM. These are determined based on the amplitude of change in growth rate of species in the presence and absence of another species in the community (>10% of change in growth of the particular species when in the presence of another species relative to the absence of another species indicates significant interaction), and the sign of the change (positive or negative). The information about the calculations of change and the type of interaction predicted in each community is added to the original table with the growth rates.
    :param inGRs: path to the file with the table listing the growth rates of the two species in a two-species community metabolic model in the presence and absence of another species in the community.
    :param outInter: path to the file that will contain the information contained in the file with growth rates, plus information regarding the the types of interactions predicted to be occurring in the community
    :return outInter: file with the interactions that are predicted to be occurring between species in a two-species community.
    '''

    cherrypy.log("We will use the information on the growth rates of the species in file %s to determine what kind of interaction is occurring between the organisms. We will output the table of interactions to %s. We will also count how many instances of each type of interaction are found" %(inGRs,outInter))

    grFile = open(inGRs, 'r')

    interactionsTableFile = open(outInter,'w')
    
    print>> interactionsTableFile, 'Model','\t','GenomeIDSpeciesA', '\t','GenomeIDSpeciesB','\t','GRSpeciesAFull','\t','GRSpeciesBFull','\t','GRASolo','\t','GRBSolo','\t','PercentChangeRawA','\t','PercentChangeRawB','\t', 'TypeOfInteraction'
    
    #next(grFile)

    # We will count how many times each interaction is predicted to occur. This information is shown in the terminal window and in the logError file.
    countMutualism = 0
    countParasitism = 0
    countCommensalism = 0
    countCompetition = 0
    countAmensalism = 0
    countNeutralism = 0



    for item in grFile:
        
        item = item.replace("_",".")
        item = item.replace("A.","")
        item = item.replace(".model","")


        try:
            itemNew = item.split('\t')
            item = item.rstrip()

            # Calculation of the effect of the presence of a competing species in the growht rate of species A
            if float(itemNew[5]) != 0:
                percentChangeRawA = (float(itemNew[3])-float(itemNew[5]))/float(itemNew[5])
            else:
                percentChangeRawA = (float(itemNew[3])-float(itemNew[5]))/float(1e-25)


            # Calculation of the effect of the presence of a competing species in the growht rate of species B
            if float(itemNew[6]) != 0:
                percentChangeRawB = (float(itemNew[4])-float(itemNew[6]))/float(itemNew[6])
            else:
                percentChangeRawB = (float(itemNew[4])-float(itemNew[6]))/float(1e-25)
        
            # Assign a type of interaction to the community based on the percent change
            #cherrypy.log('For this item, species A changed %f percent, and species B changed %f percent.'%(percentChangeRawA,percentChangeRawB))

            if percentChangeRawA > 0.1 and percentChangeRawB > 0.1:
                typeOfInteraction = 'Mutualism'
                countMutualism += 1

            elif percentChangeRawA > 0.1 and percentChangeRawB < -0.1:
                typeOfInteraction ='Parasitism'
                countParasitism += 1

            elif percentChangeRawA > 0.1 and percentChangeRawB > -0.1 and percentChangeRawB < 0.1:
                typeOfInteraction = 'Commensalism'
                countCommensalism += 1

            elif percentChangeRawA < -0.1 and percentChangeRawB > 0.1:
                typeOfInteraction = 'Parasitism'
                countParasitism += 1

            elif percentChangeRawA < -0.1 and percentChangeRawB < -0.1:
                typeOfInteraction = 'Competition'
                countCompetition += 1

            elif percentChangeRawA < -0.1 and percentChangeRawB > -0.1 and percentChangeRawB < 0.1:
                typeOfInteraction = 'Amensalism'
                countAmensalism += 1

            elif percentChangeRawA > -0.1 and percentChangeRawA < 0.1 and percentChangeRawB > 0.1:
                typeOfInteraction = 'Commensalism'
                countCommensalism += 1

            elif percentChangeRawA > -0.1 and percentChangeRawA < 0.1 and percentChangeRawB < -0.1:
                typeOfInteraction = 'Amensalism'
                countAmensalism += 1

            elif percentChangeRawA > -0.1 and percentChangeRawA < 0.1 and percentChangeRawB > -0.1 and percentChangeRawB < 0.1:
                typeOfInteraction = 'Neutralism'
                countNeutralism += 1

            else:
                typeOfInteraction = 'Empty'
                #cherrypy.log('Attention! For model %s , an interaction was identified as Empty. Something is wrong!' %item[0])#todo print out to user

        except Exception as e:
            print e

        

        # Create the interactions table. See what the growth rates file looks like, and get the appropriate columns from there, and merge the information with the information in the file with the growth rates.

        print>>interactionsTableFile, item,'\t', percentChangeRawA,'\t', percentChangeRawB, '\t', typeOfInteraction


    # Report the counts for each interaction type.
    cherrypy.log("We finished creating the interactions table, and saved it to the file %s ." %interactionsTableFile)
    cherrypy.log("We counted %d interactions that were identified as Mutualism." %countMutualism)
    cherrypy.log("We counted %d interactions that were identified as Parasitism." %countParasitism)
    cherrypy.log("We counted %d interactions that were identified as Commensalism." %countCommensalism)
    cherrypy.log("We counted %d interactions that were identified as Competition." %countCompetition)
    cherrypy.log("We counted %d interactions that were identified as Amensalism." %countAmensalism)
    cherrypy.log("We counted %d interactions that were identified as Neutralism." %countNeutralism)

        
    interactionsTableFile.close()


# In[1]:


def apply_diet(diet, models_dir, models_dieted):
    '''
    In this function we use cobrapy to calculate the growth rates of the two species that make up the two species community metabolic models under particular metabolite availability conditions. We start by loading the community model (full model) into 3 distinct Model objects with cobrapy. We then change the fluxes of the exchange reactions of the external model so they have lower bounds corresponding to whichever 'Diet' condition the user specifies. We then run a flux balance analysis on the full model, optimizing the biomass reactions of the two species that make up the community at the same time. It then remove all reactions whose IDs start with modelA from the model in modelMinusA, thus leaving only the reactions from modelB and from the external compartment. It then runs a FBA on it, maximizing the biomass reaction for the model with the tag modelB. It then does the samething but for reactions tagged with modelB on modelMinusB. The optimal flux values for the biomass reactions of each species resulting from optimization in the full model and in each model containing only one species, which correspond to predicted growth rates, are then exported to a table in the tab-delimited text formal to a folder chosen by the user.
    :param diet: the metabolite availability conditions. Default on MMinte is complete, but the user can choose another value ('Variant1 through 10')
    :param comFolder: path to the folder containing all the two-species community metabolic models.
    :return outputGRs: table with growth rate information for each of the species belonging to a two-species community metabolic model in the presence and absence of another species.
    '''
    
    allModels = getListOfModels(models_dir)


    for item in range(len(allModels)):
        
        '''
        @summary: load the model, all versions that will be manipulated in the analysis.
        '''

        # Import the models with cobrapy
        #cherrypy.log(allModels[item])
        modelFull = cobra.io.read_sbml_model(allModels[item])


        #cherrypy.log('We successfully loaded the file %s into a Model object with id, %s. They should all have the same id.'%(allModels[item],modelFull.id))
        
        # Determine what the objective function is. It should be composed of two reactions, the biomass reactions for each of the species that compose the model. Store the biomass reactions in a new variable to be used later.
        modelID = modelFull.id
 
        
        # Open the metabolite conditions file ('Diet')
        dietValues = open(diet,'r')
    
        # Default 'Diet' for MMinte is Complete. After choosing the 'Diet' change the lower bounds for the exchange reactions of the external compartment.

        #cherrypy.log('The diet chosen for this particular run of the functions was %s' %diet)
        #cherrypy.log(str(modelFull))
        for line in dietValues:
            try:
                new_line = line.rstrip('\n').split('\t')
                modelFull.reactions.get_by_id(new_line[0][0:-3]).lower_bound = -float(new_line[1])
                cherrypy.log(new_line[0][0:-3]+"found")
                cherrypy.log("eureka")
            except:
                cherrypy.log(new_line[0][0:-3])
                continue


        #cherrypy.log('We finished changing the lower bounds for the fluxes of the exchange reactions in the models to better fit the availability of metabolites for the microbial communities we are simulating the growth of. ')
        dietValues.close()

        # Run FBA on Full model
        cobra.io.write_sbml_model(modelFull,models_dieted+modelID+".xml")




