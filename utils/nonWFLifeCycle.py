#SLiM-Tree's non-Wright-Fisher life cycle and what it implies for calibration.
#
#Life cycle (written by writeSLiM.write_reproduction and the non-WF early() events), with census size K:
#  1. Reproduction: every individual is the first parent of b = OFFSPRING_PER_TICK offspring; each offspring's second
#     parent is drawn uniformly from the subpopulation.
#  2. Survival: parents and offspring alike survive with probability w * K / sum(w), where w is fitness. The expected
#     number of survivors is K, and each individual survives with probability about 1 / (1 + b).
#
#Assumptions behind the formulas below: every individual faces survival every tick, with no age effects; there is no
#fecundity selection; mates are drawn at random; selection is genic (multiplicative per copy) and weak; K is large
#and roughly constant. Under these, with phi = b / (1 + b) the fraction of survivors that are newborns:
#  * generation time = 1 / phi = (1 + b) / b ticks (a surviving lineage passes through a birth with probability phi)
#  * two gene copies coalesce in one tick with probability [2 phi (1 - phi) + phi^2 (4b - 1) / (4b)] / (2K)
#    (adult-newborn pairs 1/(2K); newborn pairs 1/(2K) unless both copies come from their first parents, in which
#    case they coalesce only if the newborns are siblings; adult pairs never), so the effective size per tick is
#        Ne_tick = 4K (1 + b)^2 / (b (4b + 7))
#    and per generation Ne_gen = phi * Ne_tick = 4K (1 + b) / (4b + 7), giving neutral diversity 4 Ne_gen mu
#  * survival selection changes an allele's frequency by s p (1 - p) per tick, so scaled selection is 4 Ne_tick s
#
#The fitness profiles are fitted for a Wright-Fisher population of K (stationary distribution proportional to
#exp(4 K f)), so selection coefficients are multiplied by K / Ne_tick = selection_scale(b) to restore 4 K f.
#For b = 1: generation time 2 ticks, Ne_tick = 16K/11, Ne_gen = 8K/11, selection scale 11/16. As b grows the
#life cycle approaches Wright-Fisher and all three tend to 1 (scale) or K (sizes).
#If the life cycle changes in any other way, these formulas must be re-derived.

#Offspring per individual per tick (b)
OFFSPRING_PER_TICK = 1


def generation_time(b = OFFSPRING_PER_TICK):
    return (1.0 + b) / b


def ne_per_tick_ratio(b = OFFSPRING_PER_TICK):
    #Ne per tick / K
    return 4.0 * (1 + b) ** 2 / (b * (4 * b + 7))


def ne_per_generation_ratio(b = OFFSPRING_PER_TICK):
    #Ne per generation / K
    return 4.0 * (1 + b) / (4 * b + 7)


def selection_scale(b = OFFSPRING_PER_TICK):
    #K / Ne per tick: multiplies selection coefficients so that scaled selection is 4 K s, as in the fitting
    return b * (4 * b + 7) / (4.0 * (1 + b) ** 2)
