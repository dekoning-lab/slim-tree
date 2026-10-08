import re
import sys


#Converts branch lengths from substitutions per site to the simulation's time units: generations for Wright-Fisher models
#(ticks_per_generation = 1) or ticks for non-WF models, which pass their generation time in ticks and their effective
#population size per generation (see utils/nonWFLifeCycle.py and readInput.process_input)
def substitutions_to_generations(tree_str, population_size, mutation_rate, ticks_per_generation = 1.0):
    tree_str = tree_str.replace(" ", "")

    # Sum only numbers that contain a decimal point (branch lengths, not bootstrap integers)
    numbers = [float(x) for x in re.findall(r'\d+\.\d+', tree_str)]
    if not numbers:
        print("Could not find any decimal branch lengths in the tree. "
              "Ensure the tree has branch lengths in substitutions. Exiting.", file=sys.stderr)
        sys.exit(1)
    total = sum(numbers)

    # First pass: normalize so all decimal branch lengths collectively sum to 30
    def scale_to_30(match):
        return str(float(match.group()) * (30 / total))
    tree_str = re.sub(r'\d+\.\d+', scale_to_30, tree_str)

    # Compute the generation-time correction factor
    N, v = population_size, mutation_rate
    Tbfix = ((1 + (4 * N * v)) - 10.32 * N * N * v * v) / v

    # Second pass: scale to generations and round to integers
    def to_generations(match):
        return str(round(float(match.group()) * (Tbfix / 3) * ticks_per_generation))
    return re.sub(r'\d+\.\d+', to_generations, tree_str)
