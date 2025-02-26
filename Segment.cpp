#include "Segment.h"



std::vector<std::tuple<uint64_t, uint64_t, double, std::string>> SegmentDNACostAndWord(
	const std::string& sequence,
	int minSegmentSize,
	int wordSize,
	int lookaheadSize/*,
	const std::function<std::pair<double, std::string>(std::string_view, int)>& costFunction*/
)
{
	if (sequence.size() < static_cast<size_t>(minSegmentSize * wordSize)) {
		throw std::invalid_argument("Sequence length must be at least the minimum segment size in words.");
	}

	std::vector<std::tuple<uint64_t, uint64_t, double, std::string>> segments;
	uint64_t currentStart = 0; // Starting position of the current segment
	uint64_t n = sequence.size();
	std::vector<std::vector<int>> leftMatrix, rightMatrix;



	while (currentStart < n)
	{
		if (!rightMatrix.empty())
		{
			leftMatrix = rightMatrix;
			rightMatrix.clear();
		}
		double bestScore = -1.0; // Track the best score in the current lookahead window
		uint64_t bestEnd = -1;        // Track the ending position of the best segment
		std::pair<double, std::string> bestSegment;
		int bestLookaheadStep = 0; // Track how far into the lookahead the max was found

		// Initialize the left and right segment sizes
		int leftSegmentSizeTemp = minSegmentSize * wordSize;
		uint64_t leftSegmentSize = leftSegmentSizeTemp;
		uint64_t rightSegmentSize = minSegmentSize * wordSize;

		// Loop through the lookahead steps
		for (int i = 0; i < lookaheadSize; ++i)
		{
			// Define the current left and right segments
			std::string_view leftSegment(sequence.data() + currentStart, leftSegmentSize);

			uint64_t currentEnd = currentStart + leftSegmentSize;
			if (currentEnd + rightSegmentSize > n) break; // Prevent out-of-bound errors
			std::string_view rightSegment(sequence.data() + currentEnd, rightSegmentSize);

			if (leftMatrix.empty() )
			{
				leftMatrix = GenerateOccurrenceMatrix(leftSegment, wordSize);
			}
			if (rightMatrix.empty())
			{
				rightMatrix = GenerateOccurrenceMatrix(rightSegment, wordSize);
			}
			if (i != 0)
			{
				int segmentSize = ((wordSize * 2) - 1);
				//Left Matrix Calculation
				uint64_t startOfSegmentToAdd = (currentStart + leftSegmentSize) - segmentSize;
					std::string_view leftSegmentToAdd(sequence.data() + startOfSegmentToAdd, segmentSize);
					auto leftMatrixToAdd = GenerateOccurrenceMatrix(leftSegmentToAdd, wordSize);
					leftMatrix = sumMatrices(leftMatrix, leftMatrixToAdd);
				//COST FUNC HERE
					uint64_t startOfSegmentToRemove = (currentEnd - wordSize);

				//Right Matrix Calculation
					std::string_view leftSegmentToRemove(sequence.data() + startOfSegmentToRemove, segmentSize);
					auto leftMatrixToRemove = GenerateOccurrenceMatrix(leftSegmentToRemove, wordSize);
					rightMatrix = subtractMatrices(rightMatrix, leftMatrixToRemove);

					uint64_t startOfSegmentToAddFromRight = (currentEnd + rightSegmentSize) - segmentSize;
					std::string_view rightSegmentToAdd(sequence.data() + startOfSegmentToAddFromRight, segmentSize);
					auto rightMatrixToAdd = GenerateOccurrenceMatrix(rightSegmentToAdd, wordSize);
					rightMatrix = sumMatrices(rightMatrix, rightMatrixToAdd);
					//Cost Function HERE
			}
			

			// Calculate costs for the left and right segments
			auto left = CalculatePercentageSumAndWord(leftMatrix);
			auto right = CalculatePercentageSumAndWord(rightMatrix);

			double totalScore = left.first + right.first;

			// Check if this score is the best so far
			if (totalScore > bestScore) {
				bestScore = totalScore;
				bestEnd = currentEnd;
				bestSegment = left;
				bestLookaheadStep = i;
			}

			// Adjust the left and right segment sizes for the next iteration
			leftSegmentSize += wordSize; // Add one word to the left segment
			currentEnd += wordSize;      // Move the right segment one word forward
		}

		// Backtrace if the maximum was found before the end of the lookahead loop
		if (bestLookaheadStep < lookaheadSize - 1) {
			// Backtracking logic can be expanded as needed
			// For now, simply record the best segment and move on
		}

		// Add the best left segment to the results and move the current start position
		if (bestEnd != -1) {
			segments.emplace_back(currentStart, bestEnd, bestSegment.first, bestSegment.second);
			currentStart = bestEnd; // Move the start position to the end of the best segment
		}
		else {
			break; // No valid segments found, terminate
		}
	}

	return segments;
}