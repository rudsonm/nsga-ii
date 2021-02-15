#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <tgmath.h>
#include <map>
#include <cmath>
#include <iostream>

#define PENALTY_FACTOR 10

struct Solution {
    int id = 0;
    std::vector<int> requirementTeam;
    std::vector<int> requirementSequence;
    std::vector<float> objectives;
    std::vector<Solution*> dominates;
    float penalty = 1.;
    float distance = 0.;
    int rank;

    Solution() {}

    Solution(int solutionId, int numRequirements, int numTeams) {
        id = solutionId;
        requirementTeam = std::vector<int>(numRequirements, 0);
        requirementSequence = std::vector<int>(numTeams, 0);
    }

    float getCost() {
        return objectives.at(0) + PENALTY_FACTOR * penalty;
    }

    float getSatisfaction() {
        return objectives.at(1) - PENALTY_FACTOR * penalty;
    }
};

struct Instance {
    std::vector<float> customerWeight;
    std::vector<std::vector<float>> requirementImportance;
    std::vector<std::vector<float>> requirementCost;
    std::vector<float> teamHourCapacity;
    std::vector<float> teamHourCost;
    std::vector<std::pair<int, int>> prerequisites;

    void print() {
        printf("-------------INSTANCE DETAILS-------------\n");

        printf("CUSTOMERS\n");
        for(int i = 0; i < customerWeight.size(); i++) {
            printf("C: %i\tW: %2.f\n", i, customerWeight.at(i));
        }

        printf("TEAMS\n");
        for(int i = 0; i < teamHourCapacity.size(); i++) {
            printf("T: %i\tCap: %2.fh\tCost: $%2.f per hour\n", i, teamHourCapacity.at(i), teamHourCost.at(i));
        }

        printf("REQUIREMENTS\n");        
        printf("Cost in hours per team\n");
        for(int i = 0; i < teamHourCapacity.size(); i++) {
            printf("\tT%i", i);
        }
        for(int i = 0; i < requirementCost.size(); i++) {
            printf("\nR%i", i);
            for(int j = 0; j < requirementCost.at(i).size(); j++) {
                printf("\t%2.fh", requirementCost.at(i).at(j));
            }
        }
        
        printf("\nImportance per client\n");
        for(int i = 0; i < customerWeight.size(); i++) {
            printf("\tC%i", i);
        }
        for(int i = 0; i < requirementImportance.size(); i++) {
            printf("\nR%i", i);
            for(int j = 0; j < requirementImportance.at(i).size(); j++) {
                printf("\t%2.f", requirementImportance.at(i).at(j));
            }
        }
    }
};

struct ConstraintValidator {
    int getTotalInfractions(const Instance &instance, const Solution &solution) {
        return teamCapacity(instance, solution) +
               teamSequence(instance, solution);
    }

    int teamCapacity(const Instance &instance, const Solution &solution) {
        std::map<int, float> teamAlloc;
        for (int requirement = 0; requirement < solution.requirementTeam.size(); requirement++) {
            int team = solution.requirementTeam.at(requirement);
            if (team == 0)
                continue;
            float cost = instance.requirementCost.at(requirement).at(team-1);
            if (teamAlloc.find(team) == teamAlloc.end()) {
                teamAlloc.insert(std::make_pair(team, cost));
            } else {
                teamAlloc.at(team) += cost;
            }
        }
        int infractions = 0;
        for (auto const& x : teamAlloc) {
            int team = x.first;
            int cost = x.second;
            if (cost > instance.teamHourCapacity.at(team-1))
                infractions++;
        }
        return infractions;
    }

    int teamSequence(const Instance &instance, const Solution &solution) {
        int infractions = 0;
        std::map<int, std::vector<int>> reqByTeam;
        for (int requirement = 0; requirement < solution.requirementTeam.size(); requirement++) {
            int team = solution.requirementTeam.at(requirement);
            if (team == 0)
                continue;
                
            if (reqByTeam.find(team) == reqByTeam.end())
                reqByTeam.insert(std::make_pair(team, std::vector<int>({requirement})));
            else
                reqByTeam.at(team).push_back(requirement);
        }        
        for (auto const &prerequisite : instance.prerequisites) {
            int reqI = prerequisite.first,
                reqJ = prerequisite.second;
            
            int teamI = solution.requirementTeam.at(reqI);
            float startI = 0.;
            for (int requirement : reqByTeam.at(teamI))
                if (solution.requirementSequence.at(requirement) < solution.requirementSequence.at(reqI))
                    startI += instance.requirementCost.at(requirement).at(teamI-1);
            
            float startJ = 0.;
            float teamJ = solution.requirementTeam.at(reqJ);
            for (int requirement : reqByTeam.at(teamJ))
                if (solution.requirementSequence.at(requirement) < solution.requirementSequence.at(reqJ))
                    startI += instance.requirementCost.at(requirement).at(teamJ-1);
            
            float durationI = instance.requirementCost.at(reqI).at(teamJ-1);
            if (startI + durationI > startJ)
                infractions++;
        }
        return infractions;
    }
};

Instance generateRandomInstance(const int &numCustomer, const int &numRequirement, const int &numTeam) {
    Instance instance;
    for (int i = 0; i < numCustomer; i++) {
        instance.customerWeight.push_back(rand() % 3 + 1);
    }

    for (int i = 0; i < numTeam; i++) {
        instance.teamHourCapacity.push_back(rand() % 20 + 20);
        instance.teamHourCost.push_back(rand() % 15 + 20);
    }

    for (int i = 0; i < numRequirement; i++) {
        instance.requirementImportance.push_back(std::vector<float>(numCustomer, 0.));
        for (int j = 0; j < numCustomer; j++) {
            instance.requirementImportance.at(i).at(j) = rand() % 10 + 1;
        }
        
        instance.requirementCost.push_back(std::vector<float>(numTeam, 0.));
        for (int j = 0; j < numTeam; j++) {
            instance.requirementCost.at(i).at(j) = (1 - instance.teamHourCost.at(j) / 35) * (rand() % 20 + 1);
        }
    }

    return instance;
}

Solution generateRandomSolution(const int &id, const int &numRequirements, const int &numTeams) {
    Solution solution(id, numRequirements, numTeams);

    for (int &sequence : solution.requirementSequence)
        sequence = std::rand() % (numRequirements + 1);
    for (int &team : solution.requirementTeam)
        team = std::rand() % (numTeams + 1);

    return solution;
}

std::vector<Solution> generateRandomPopulation(const int &numPopulation, const int &numRequirements, const int &numTeams) {
    std::vector<Solution> solutions;

    for (int i = 0; i < numPopulation; i++) {
        Solution randomSolution = generateRandomSolution(i, numRequirements, numTeams);
        solutions.push_back(randomSolution);
    }
    return solutions;
}

void assignCostValue(Instance instance, Solution &solution) {
    float costValue = 0.;
    for (int requirement = 0; requirement < solution.requirementTeam.size(); requirement++) {
        int team = solution.requirementTeam[requirement] ;
        if (team == 0)
            continue;
        float teamCost = instance.teamHourCost.at(team - 1);
        float requirementCost = instance.requirementCost.at(requirement).at(team - 1);
        costValue += requirementCost * teamCost;
    }
    solution.objectives.push_back(costValue);
}

void assignSatisfactionValue(Instance instance, Solution &solution) {
    float satisfactionValue = 0.;
    for (int requirement = 0; requirement < solution.requirementTeam.size(); requirement++) {
        int team = solution.requirementTeam[requirement];
        if (team == 0)
            continue;
        for (int client = 0; client < instance.requirementImportance.at(requirement).size(); client++) {
            float satisfaction = instance.requirementImportance.at(requirement).at(client);
            float weight = instance.customerWeight.at(client);
            satisfactionValue += satisfaction * weight;
        }        
    }
    solution.objectives.push_back(satisfactionValue);
}

void evaluateSolutions(Instance instance, std::vector<Solution> &solutions) {
    std::vector<std::pair<int, float>> costValues;
    std::vector<std::pair<int, float>> satisfactionValues;
    ConstraintValidator cv;
    for (int i = 0; i < solutions.size(); i++) {
        if (solutions.at(i).objectives.size() > 0)
            solutions.at(i).objectives.clear();
        assignCostValue(instance, solutions.at(i));
        assignSatisfactionValue(instance, solutions.at(i));
        solutions.at(i).penalty = cv.getTotalInfractions(instance, solutions.at(i));
    }
}

bool maxCompare(std::pair<int, float> a, std::pair<int, float> b) {
    return a.second < b.second;
}

bool minCompare(std::pair<int, float> a, std::pair<int, float> b) {
    return b.second < a.second;
}

std::vector<Solution> binaryTournamentSelection(int n, std::vector<Solution> parents) {
    std::vector<Solution> winners;
    int randomIndexA, randomIndexB;
    do {
        do {
            randomIndexA = rand() % parents.size();
            randomIndexB = rand() % parents.size();
        } while (randomIndexA != randomIndexB);
        
        Solution parentA = parents.at(randomIndexA),
                 parentB = parents.at(randomIndexB);
        
        Solution winner;
        int winnerIndex;
        if (parentA.rank < parentB.rank) {
            winner = parentA;
            winnerIndex = randomIndexA;
        } else if (parentA.rank == parentB.rank) {
            winner = parentA.distance > parentB.distance ? parentA : parentB;
            winnerIndex = parentA.distance > parentB.distance ? randomIndexA : randomIndexB;
        } else {
            winner = parentB;
            winnerIndex = randomIndexB;
        }
        winners.push_back(winner);
        parents.erase(parents.begin() + winnerIndex);
    } while(winners.size() < n && parents.size() > 1);
    return winners;
}

std::vector<Solution> kPointCrossover(int n, int k, std::vector<Solution> parents, int numRequirement, int numTeam) {
    std::vector<Solution> childrens;
    do {
        int randomIndexA, randomIndexB;
        do {
            randomIndexA = rand() % parents.size();
            randomIndexB = rand() % parents.size();
        } while (randomIndexA != randomIndexB);
        
        Solution parentA = parents.at(randomIndexA),
                 parentB = parents.at(randomIndexB);
        
        Solution children1, children2;
        int batch = parentA.requirementSequence.size() / k;
        for (int i = 0; i < batch; i++) {
            Solution parent1 = i % 2 == 0 ? parentA : parentB,
                     parent2 = i % 2 == 0 ? parentB : parentA;
            int start = batch * i,
                end = i == batch - 1 ? parentA.requirementSequence.size() : start + batch;
            for (int j = start; j < end; j++) {
                children1.requirementSequence.push_back(parent1.requirementSequence.at(j));
                children1.requirementTeam.push_back(parent1.requirementTeam.at(j));

                children2.requirementSequence.push_back(parent2.requirementSequence.at(j));
                children2.requirementTeam.push_back(parent2.requirementTeam.at(j));
            }
        }
        childrens.push_back(children1);
        childrens.push_back(children2);
    } while(childrens.size() < n);
    return childrens;
}

std::vector<std::vector<int>> getParetoFrontsAndAssignSolutionRank(std::vector<Solution> &solutions) {
    std::vector<int> dominationCount(solutions.size(), 0);
    std::vector<std::vector<int>> fronts(1, std::vector<int>());    
    for (int i = 0; i < solutions.size(); i++) {
        for (int j = 0; j < solutions.size(); j++) {
            if (i == j)
                continue;
            bool iDominatesJ = solutions.at(i).getCost() < solutions.at(j).getCost()
                            && solutions.at(i).getSatisfaction() > solutions.at(j).getSatisfaction();
            bool jDominatesI = solutions.at(j).getCost() < solutions.at(i).getCost()
                            && solutions.at(j).getSatisfaction() > solutions.at(i).getSatisfaction();
            if (iDominatesJ)
                solutions.at(i).dominates.push_back(&solutions.at(j));
            else if(jDominatesI)
                dominationCount.at(i)++;
        }
        if (dominationCount.at(i) == 0) {
            solutions.at(i).rank = 0;
            fronts.at(0).push_back(i);
        }
    }
    while(fronts.back().size() != 0) {
        std::vector<int> nextFront;
        for (int i = 0; i < fronts.back().size(); i++) {
            Solution currentSolution = solutions.at(fronts.back().at(i));
            for (int j = 0; j < currentSolution.dominates.size(); j++) {
                int jIndex = currentSolution.dominates.at(j)->id;
                dominationCount.at(jIndex)--;
                if (dominationCount.at(jIndex) == 0) {
                    solutions.at(jIndex).rank = fronts.size();
                    nextFront.push_back(jIndex);
                }
            }
        }
        fronts.push_back(nextFront);
    }
    return fronts;
}

std::vector<float> assignCrowdingDistance(std::vector<Solution> &solutions) {
    std::vector<std::pair<int, float>> costValues;
    std::vector<std::pair<int, float>> satisfactionValues;
    printf("A");
    for (int i = 0; i < solutions.size(); i++) {
        std::pair<int, float> pairValue = std::make_pair(i, solutions.at(i).getCost());
        costValues.push_back(pairValue);
        pairValue = std::make_pair(i, solutions.at(i).getSatisfaction());
        satisfactionValues.push_back(pairValue);
    }
    printf("B");
    std::sort(costValues.begin(), costValues.end(), minCompare);
    std::pair<float,float>costValuesMinMax(costValues.front().second, costValues.back().second);    
    
    std::sort(satisfactionValues.begin(), satisfactionValues.end(), maxCompare);
    std::pair<float, float>satisfactionValuesMinMax(satisfactionValues.front().second, satisfactionValues.back().second);
    printf("C");
    for (int i = 1; i < solutions.size() - 1; i++) {
        int solutionIndex = costValues.at(i).first;
        float distance = (costValues.at(i + 1).second - costValues.at(i - 1).second) / (costValuesMinMax.first - costValuesMinMax.second);
        solutions.at(solutionIndex).distance += std::fabs(distance);

        solutionIndex = satisfactionValues.at(i).first;
        distance = (satisfactionValues.at(i + 1).second - satisfactionValues.at(i - 1).second) / (satisfactionValuesMinMax.first - satisfactionValuesMinMax.second);
        solutions.at(solutionIndex).distance += std::fabs(distance);
    }
    printf("D");
}

std::vector<Solution> mergeParentsWithChilds(std::vector<Solution> parents, std::vector<Solution> childs) {
    std::vector<Solution> population;
    int id = 0;
    for (Solution parent : parents) {
        parent.id = id++;
        parent.rank = 0;
        parent.distance = 0.;
        parent.dominates.clear();
        population.push_back(parent);
    }
    for (Solution child : childs) {
        child.id = id++;
        population.push_back(child);
    }
    return population;
}

int main() {
    const int NUM_GENERATIONS = 5,
              POPULATION_SIZE = 25,
              NUM_CUSTOMERS = 5,
              NUM_REQUIREMENTS = 8,
              NUM_TEAMS = 3;

    Instance instance = generateRandomInstance(NUM_CUSTOMERS, NUM_REQUIREMENTS, NUM_TEAMS);    
    instance.print();

    std::vector<Solution> solutions = generateRandomPopulation(POPULATION_SIZE, NUM_REQUIREMENTS, NUM_TEAMS);
    for (int g = 0; g < NUM_GENERATIONS; g++) {
        std::cout << std:: endl <<  "GEN: " << g + 1 << std::endl;
        evaluateSolutions(instance, solutions);
        std::vector<std::vector<int>> fronts = getParetoFrontsAndAssignSolutionRank(solutions);
        for (int i = 0; i < fronts.size(); i++)
            for (int j = 0; j < fronts.at(i).size(); j++)
                std::cout << "Front: " << i << ", Sol: " << solutions.at(fronts.at(i).at(j)).id << std::endl;
        std::vector<Solution> nextGen;
        int i = 0;
        while(nextGen.size() < POPULATION_SIZE && i < fronts.size()) {
            std::vector<Solution> frontSolutions;
            for (int j = 0; j < fronts.at(i).size(); j++) {
                int index = fronts.at(i).at(j);
                frontSolutions.push_back(solutions.at(index));
            }
            int n = fronts.at(i).size();
            if (fronts.at(i).size() + nextGen.size() > POPULATION_SIZE) {
                n = POPULATION_SIZE - nextGen.size();
                assignCrowdingDistance(frontSolutions);
            }
            nextGen.insert(nextGen.begin(), frontSolutions.begin(), frontSolutions.begin() + n);
            i++;
        }
        std::vector<Solution> parents = binaryTournamentSelection(POPULATION_SIZE / 2, nextGen);
        std::vector<Solution> childs = kPointCrossover(POPULATION_SIZE, 4, parents, NUM_REQUIREMENTS, NUM_TEAMS);
        solutions = mergeParentsWithChilds(parents, childs);
    }

    return 0;
}