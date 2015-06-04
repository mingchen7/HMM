
# Function: f(x) = abs(x^2 + y^2 + z^2)
# x, y, z = [-10, 10]

EvaluationFunction = function(x, y, z)
{
    return(abs(x^2 + 2 * y^2 + 3 * z^2));
}


# example: GA.MainFunction();
GA.MainFunction = function()
{
        InitialGen = matrix(c(0,1,1,
                          -1,-2,-2,
                          -2,-1,1,
                          -3,-5,6,
                          -1,4,2,
                          -2,7,1,
                          3,-2,-1,
                          5,-4,-1,
                          4,6,-2,
                          3,2,-1,
                          1,2,3,
                          1,2,4,
                          -1,-2,-3,
                          -1,-2,-4,
                          -2,-2,-1,
                          -2,-2,-3,
                          -3,-1,-2,
                          3,2,1,
                          2,2,2,
                          2,1,2), 
                        ncol = 3, byrow = TRUE);
    
    for(i in 1:2000)
    {
        InitialGen = GA.NewGeneration(InitialGen);
        print(EvaluationFunction(InitialGen[1, 1], InitialGen[1, 2], InitialGen[1, 3]));
    }
    
    print(InitialGen[1, ]);
}



GA.NewGeneration = function(originalGen)
{
    NewGeneration = c(); #matrix(data = -999, nrow = nrow(originalGens), ncol = ncol(originalGens));
    
    SortedFitness = GA.Evaluate(originalGen);
    SortedOriginalGen = GA.SortGen(originalGen,SortedFitness$ix);
    
    NewGeneration = GA.Elitism(SortedOriginalGen);
    
    for(i in 1:((nrow(originalGen))/ 2 - 1))
    {
        SelectedChromos = GA.SelectedChromosomes(SortedOriginalGen);
        CrossedChromos = GA.Crossover(SelectedChromos);
        MutatedChromos = GA.Mutation(CrossedChromos);
        
        NewGeneration = rbind(NewGeneration, MutatedChromos);
    }
    
    return(NewGeneration);
}

GA.Evaluate = function(originalGen)
{
    Fitness = c();
    
    for(i in 1: nrow(originalGen))
    {
        Fitness[i] = EvaluationFunction(originalGen[i, 1], originalGen[i, 2], originalGen[i, 3]);
    }
    
    SortedFitness = sort(Fitness, index.return = TRUE);
    
    return(SortedFitness);
}


GA.Elitism = function(sortedGeneration) # from min to max
{
    return(sortedGeneration[c(1,2), ]);
}

GA.SortGen = function(originalGen, sortIndex)
{
    return(originalGen[sortIndex, ]);
}

GA.SelectedChromosomes = function(sortedGeneration)
{
    FirstChromoIndex = GA.RankSelection(sortedGeneration);
    SecondChromoIndex = GA.RankSelection(sortedGeneration);
    
    #print(sprintf('%02d, %02d', FirstChromoIndex, SecondChromoIndex));
    
    return(sortedGeneration[c(FirstChromoIndex, SecondChromoIndex), ]);
}

GA.Mutation = function(crossedChromosomes)
{
    MutationProbability = 0.03;
    RandomProbability = runif(1);
    
    if(RandomProbability <= MutationProbability)
    {
        RandomVar = runif(1);
        if(RandomVar <= 1/3) # change  [x,y,z], or [x, y] or [x]
        {
            for(j in 1: sample(1:3, 1))
            {
                crossedChromosomes[2, j] = crossedChromosomes[2, j] + Generate.Sign() * 0.01;
            }
        }
        else if(RandomVar > 2/3) # change  [y,z], or [y]
        {
            for(j in 2: sample(2:3, 1))
            {
                crossedChromosomes[2, j] = crossedChromosomes[2, j] + Generate.Sign() * 0.01;
            }
        }
        else # change [z]
        {
            crossedChromosomes[2, 3] = crossedChromosomes[2, 3] + Generate.Sign() * 0.01;
        }
    }
    
    return(crossedChromosomes); 
    
}

Generate.Sign = function()
{
    RandomVar = runif(1, min = -1, max = 1);
    
    return(RandomVar / abs(RandomVar));
}

GA.Crossover = function(selectedChromosomes)
{
    CrossoverProbability = 0.8;
    RandomProbability = runif(1);
    
    if(RandomProbability <= CrossoverProbability)
    {
        #Two types of crossover: [1]start from y, switch y and z; [2]start from z, switch z only
        Option = runif(1);
        if(Option <= 0.5) #[1]start from y, switch y and z;
        {
            TempY = selectedChromosomes[1, 2];
            TempZ = selectedChromosomes[1, 3];
            
            selectedChromosomes[1, 2] = selectedChromosomes[2, 2];
            selectedChromosomes[1, 3] = selectedChromosomes[2, 3];
            
            selectedChromosomes[2, 2] = TempY;
            selectedChromosomes[2, 3] = TempZ;
        }
        else
        {
            TempZ = selectedChromosomes[1, 3];

            selectedChromosomes[1, 3] = selectedChromosomes[2, 3];

            selectedChromosomes[2, 3] = TempZ;
        }
    }

    return(selectedChromosomes);
    
}

# example: GA.RankSelection(InitialGen)
GA.RankSelection = function(generation)
{
    LabelVector = c();
    for(i in 1: nrow(generation))
    {
        for(j in (sum(1: (i-1)) + 1): sum(1:i))
        {
            LabelVector[j] = i;
        }
    }
    
    RandomVar = runif(1);
    
    By = 1 / sum(1: nrow(generation));
    Index = ceiling(RandomVar / By);
    
    return(nrow(generation) - LabelVector[Index] + 1);
}









