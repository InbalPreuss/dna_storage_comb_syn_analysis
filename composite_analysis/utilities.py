import math
import Levenshtein


def is_within_levenshtein_distance(seq: str, target_seq: str, max_distance=4) -> tuple[int, int, int]:
    distance, start_pos_temp, end_pos = math.inf, math.inf, math.inf
    seq_len = len(seq)
    target_seq_len = len(target_seq)

    for start_pos in range(target_seq_len - seq_len + 1):
        subseq = target_seq[start_pos:start_pos + seq_len]
        distance = Levenshtein.distance(seq, subseq)
        if distance <= max_distance:
            max_distance = distance
            start_pos_temp = start_pos

    return max_distance, start_pos_temp, start_pos_temp + seq_len
