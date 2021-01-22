#include <cstdlib>
#include <vector>
#include <algorithm>
#include <limits>

struct Solution {
    std::vector<int> requirementTeam;
    std::vector<int> requirementSequence;

    Solution(int numRequirements, int numTeams) {
        requirementTeam = std::vector<int>(numRequirements, 0);
        requirementSequence = std::vector<int>(numTeams, 0);
    }
};

struct Instance {
    std::vector<float> customerWeight;
    std::vector<std::vector<float>> requirementImportance;
    std::vector<std::vector<float>> requirementCost;
    std::vector<float> teamHourCost;
};

Solution generateRandomSolution(const int numRequirements, const int numTeams) {
    Solution solution(numRequirements, numTeams);

    for (int &sequence : solution.requirementSequence)
        sequence = std::rand() % (numRequirements + 1);
    for (int &team : solution.requirementTeam)
        team = std::rand() % (numTeams + 1);

    return solution;
}

float getCostValue(Instance instance, Solution solution) {
    float costValue = 0.;
    for (int requirement = 0; requirement < solution.requirementTeam.size(); requirement++) {
        int team = solution.requirementTeam[requirement];
        if (team == 0)
            continue;
        float teamCost = instance.teamHourCost.at(team);
        float requirementCost = instance.requirementCost.at(requirement).at(team);
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

int main() {
    const int POPULATION_SIZE = 10, NUM_REQUIREMENTS = 6, NUM_TEAMS = 2;

    Instance instance;

    // GENERATE RANDOM POPULATION
    std::vector<Solution> solutions;
    for (int i = 0; i < POPULATION_SIZE; i++)
        solutions.push_back(generateRandomSolution(NUM_REQUIREMENTS, NUM_TEAMS));

    // EVALUATE OBJECTIVE VALUES
    std::vector<std::pair<int, float>> costValues;
    std::vector<std::pair<int, float>> satisfactionValues;
    for (int i = 0; i < solutions.size(); i++) {    
        Solution solution = solutions.at(i);
        costValues.push_back(std::make_pair(i, getCostValue(instance, solution)));
        satisfactionValues.push_back(std::make_pair(i, getSatisfactionValue(instance, solution)));
    }

    // CROWDING DISTANCE ASSIGNMENT
    std::sort(costValues.begin(), costValues.end(), minCompare);
    std::pair<float,float>costValuesMinMax(costValues.front().second, costValues.back().second);
    
    std::sort(satisfactionValues.begin(), satisfactionValues.end(), maxCompare);
    std::pair<float, float>satisfactionValuesMinMax(satisfactionValues.front().second, satisfactionValues.back().second);

    std::vector<float> solutionsDistance(solutions.size(), 0.);
    solutionsDistance.front() = solutionsDistance.back() = std::numeric_limits<float>::max();    
    for (int i = 1; i < solutions.size() - 1; i++) {
        int solutionIndex = costValues.at(i).first;
        solutionsDistance.at(solutionIndex) += (costValues.at(i + 1).second - costValues.at(i - 1).second) / (costValuesMinMax.first - costValuesMinMax.second);

        solutionIndex = satisfactionValues.at(i).first;
        solutionsDistance.at(solutionIndex) += (satisfactionValues.at(i + 1).second - satisfactionValues.at(i - 1).second) / (satisfactionValuesMinMax.first - satisfactionValuesMinMax.second);
    }
    return 0;
}