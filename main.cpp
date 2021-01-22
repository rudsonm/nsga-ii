#include <vector>
#include <cstdlib>

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

int main() {
    const int POPULATION_SIZE = 10, NUM_REQUIREMENTS = 6, NUM_TEAMS = 2;

    Instance instance;

    // GENERATE RANDOM POPULATION
    std::vector<Solution> solutions;
    for (int i = 0; i < POPULATION_SIZE; i++)
        solutions.push_back(generateRandomSolution(NUM_REQUIREMENTS, NUM_TEAMS));

    // EVALUATE OBJECTIVE VALUES
    std::vector<float> costValues(solutions.size(), 0.);
    std::vector<float> satisfactionValues(solutions.size(), 0.);
    for (int i = 0; i < solutions.size(); i++) {
        Solution solution = solutions.at(i);
        costValues.at(i) = getCostValue(instance, solution);
        satisfactionValues.at(i) = getSatisfactionValue(instance, solution);
    }

    return 0;
}