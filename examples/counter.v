module counter(
	input clk,
	input reset,
	output [3:0] count
);
	always @ (posedge clk) begin
		if (reset)
			count <= 4'd0;
		else
			count <= count + 4'd1;
	end

endmodule
