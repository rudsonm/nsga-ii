#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <tgmath.h>

struct Solution {
    int id = 0;
    std::vector<int> requirementTeam;
    std::vector<int> requirementSequence;
    std::vector<Solution*> dominates;
    float distance = 0.;
    int rank;    

    Solution() {}

    Solution(int solutionId, int numRequirements, int numTeams) {
        id = solutionId;
        requirementTeam = std::vector<int>(numRequirements, 0);
        requirementSequence = std::vector<int>(numTeams, 0);
    }

};

struct Instance {
    std::vector<float> customerWeight;
    std::vector<std::vector<float>> requirementImportance;
    std::vector<std::vector<float>> requirementCost;
    std::vector<float> teamHourCapacity;
    std::vector<float> teamHourCost;    
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

float getCostValue(Instance instance, Solution solution) {
    float costValue = 0.;
    for (int requirement = 0; requirement < solution.requirementTeam.size(); requirement++) {
        int team = solution.requirementTeam[requirement] ;
        if (team == 0)
            continue;
        float teamCost = instance.teamHourCost.at(team - 1);
        float requirementCost = instance.requirementCost.at(requirement).at(team - 1);
        costValue += requirementCost * teamCost;
    }
    return costValue;
}

float getSatisfactionValue(Instance instance, Solution solution) {
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
    return satisfactionValue;
}

bool maxCompare(std::pair<int, float> a, std::pair<int, float> b) {
    return a.second < b.second;
}

bool minCompare(std::pair<int, float> a, std::pair<int, float> b) {
    return b.second < a.second;
}

Solution binaryTournamentSelection(const std::vector<Solution> &parents) {
    int randomIndexA, randomIndexB;
    do {
        randomIndexA = rand() % parents.size();
        randomIndexB = rand() % parents.size();
    } while (randomIndexA != randomIndexB);
    
    Solution parentA = parents.at(randomIndexA),
             parentB = parents.at(randomIndexB);
    
    if (parentA.rank < parentB.rank) {
        return parentA;
    } else if (parentA.rank == parentB.rank) {
        return parentA.distance > parentB.distance ? parentA : parentB;
    } else {
        return parentB;
    }
}

int main() {
    const int POPULATION_SIZE = 10,
              NUM_CUSTOMERS = 5,
              NUM_REQUIREMENTS = 8,
              NUM_TEAMS = 3;

    Instance instance = generateRandomInstance(NUM_CUSTOMERS, NUM_REQUIREMENTS, NUM_TEAMS);
    
    // GENERATE RANDOM POPULATION
    std::vector<Solution> solutions;
    for (int i = 0; i < POPULATION_SIZE; i++) {
        Solution randomSolution = generateRandomSolution(i, NUM_REQUIREMENTS, NUM_TEAMS);
        solutions.push_back(randomSolution);
    }
    
    // EVALUATE OBJECTIVE VALUES
    std::vector<std::pair<int, float>> costValues;
    std::vector<std::pair<int, float>> satisfactionValues;
    for (int i = 0; i < solutions.size(); i++) {    
        Solution solution = solutions.at(i);
        costValues.push_back(std::make_pair(i, getCostValue(instance, solution)));
        satisfactionValues.push_back(std::make_pair(i, getSatisfactionValue(instance, solution)));
    }

    // RANK ASSIGNMENT
    std::vector<int> dominationCount(solutions.size(), 0);
    std::vector<std::vector<Solution*>> fronts(1, std::vector<Solution*>());
    for (int i = 0; i < solutions.size(); i++) {
        for (int j = 0; j < solutions.size(); j++) {
            if (i == j)
                continue;
            bool iDominatesJ = costValues.at(i).second < costValues.at(j).second
                            && satisfactionValues.at(i).second > satisfactionValues.at(j).second;
            if (iDominatesJ) {
                solutions.at(i).dominates.push_back(&solutions.at(j));
            } else {
                dominationCount.at(i)++;
            }
        }
        if (dominationCount.at(i) == 0) {
            solutions.at(i).rank = 1;
            fronts.at(0).push_back(&solutions.at(i));
        }
    }
    while(fronts.back().size() != 0) {
        std::vector<Solution*> nextFront;
        for (int i = 0; i < fronts.back().size(); i++) {
            Solution *currentSolution = fronts.back().at(i);
            for (int j = 0; j < currentSolution->dominates.size(); j++) {
                int jIndex = currentSolution->dominates.at(j)->id;
                dominationCount.at(jIndex)--;
                if (dominationCount.at(jIndex) == 0) {
                    solutions.at(jIndex).rank++;
                    nextFront.push_back(&solutions.at(jIndex));
                }
            }
        }
        fronts.push_back(nextFront);
    }

    // CROWDING DISTANCE ASSIGNMENT
    std::sort(costValues.begin(), costValues.end(), minCompare);
    std::pair<float,float>costValuesMinMax(costValues.front().second, costValues.back().second);
    
    std::sort(satisfactionValues.begin(), satisfactionValues.end(), maxCompare);
    std::pair<float, float>satisfactionValuesMinMax(satisfactionValues.front().second, satisfactionValues.back().second);

    std::vector<float> solutionsDistance(solutions.size(), 0.);    
    for (int i = 1; i < solutions.size() - 1; i++) {
        int solutionIndex = costValues.at(i).first;
        float distance = (costValues.at(i + 1).second - costValues.at(i - 1).second) / (costValuesMinMax.first - costValuesMinMax.second);
        solutions.at(solutionIndex).distance += std::fabs(distance);

        solutionIndex = satisfactionValues.at(i).first;
        distance = (satisfactionValues.at(i + 1).second - satisfactionValues.at(i - 1).second) / (satisfactionValuesMinMax.first - satisfactionValuesMinMax.second);
        solutions.at(solutionIndex).distance += std::fabs(distance);
    }
    return 0;
}