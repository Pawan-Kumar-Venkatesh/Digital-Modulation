function run(mod_scheme,snr,num_pkt,pkt_size)

	clc;
    
    tic
	
    warning('off')

    [tx_data, tx_vec, tx_syms, rx_vec, rx_data, rx_syms, sym_err_rate, bit_err_rate, pkt_loss, tx_vec_plot, rx_vec_plot]=tx_and_rx(num_pkt, mod_scheme, pkt_size, snr);
    
    sym_err_rate
    bit_err_rate
    evm = comm.EVM();
    rmsevm=mean(evm(rx_syms,tx_syms))

    pkt_loss_percent=(pkt_loss/num_pkt)*100


    % Plotting
    plot_main(tx_vec_plot, rx_vec_plot, tx_syms, rx_syms);
    
    k=1;
    for snr=1:5:50
        % For BPSK
        mod_scheme='BPSK';
        [tx_data, tx_vec, tx_syms, rx_vec, rx_data, rx_syms, sym_err_rate, bit_err_rate, pkt_loss, tx_vec_plot, rx_vec_plot]=tx_and_rx(num_pkt, mod_scheme, pkt_size, snr);
        ber_bpsk(k)=bit_err_rate;
        ser_bpsk(k)=sym_err_rate;
        
        evm_bpsk = comm.EVM();
        rmsevm_bpsk(k)=mean(evm_bpsk(rx_syms,tx_syms));

        pkt_loss_bpsk(k)=(pkt_loss/num_pkt)*100;

        % For QPSK
        mod_scheme='QPSK';
        [tx_data, tx_vec, tx_syms, rx_vec, rx_data, rx_syms, sym_err_rate, bit_err_rate, pkt_loss, tx_vec_plot, rx_vec_plot]=tx_and_rx(num_pkt, mod_scheme, pkt_size, snr);
        ber_qpsk(k)=bit_err_rate;
        ser_qpsk(k)=sym_err_rate;
        
        evm_qpsk = comm.EVM();
        rmsevm_qpsk(k)=mean(evm_qpsk(rx_syms,tx_syms));

        pkt_loss_qpsk(k)=(pkt_loss/num_pkt)*100;

        % For 16QAM
        mod_scheme='16QAM';
        [tx_data, tx_vec, tx_syms, rx_vec, rx_data, rx_syms, sym_err_rate, bit_err_rate, pkt_loss, tx_vec_plot, rx_vec_plot]=tx_and_rx(num_pkt, mod_scheme, pkt_size, snr);
        ber_16qam(k)=bit_err_rate;
        ser_16qam(k)=sym_err_rate;

        evm_16qam = comm.EVM();
        rmsevm_16qam(k)=mean(evm_16qam(rx_syms,tx_syms));

        pkt_loss_16qam(k)=(pkt_loss/num_pkt)*100;

        % For 64QAM
        mod_scheme='64QAM';
        [tx_data, tx_vec, tx_syms, rx_vec, rx_data, rx_syms, sym_err_rate, bit_err_rate, pkt_loss, tx_vec_plot, rx_vec_plot]=tx_and_rx(num_pkt, mod_scheme, pkt_size, snr);
        ber_64qam(k)=bit_err_rate;
        ser_64qam(k)=sym_err_rate;

        evm_64qam = comm.EVM();
        rmsevm_64qam(k)=mean(evm_64qam(rx_syms,tx_syms));

        pkt_loss_64qam(k)=(pkt_loss/num_pkt)*100;

        snr_cnt(k)=snr;
        k=k+1;
    end

    figure
    plot(snr_cnt, ber_bpsk);
    hold on
    plot(snr_cnt, ber_qpsk);
    hold on
    plot(snr_cnt, ber_16qam);
    hold on
    plot(snr_cnt, ber_64qam);
    legend('BPSK','QPSK','16QAM','64QAM');
    title('Bit Error Rate');
    xlabel('SNR (dB)');
    ylabel('BER (%)');


    figure
    plot(snr_cnt, ser_bpsk);
    hold on
    plot(snr_cnt, ser_qpsk);
    hold on
    plot(snr_cnt, ser_16qam);
    hold on
    plot(snr_cnt, ser_64qam);
    legend('BPSK','QPSK','16QAM','64QAM');
    title('Symbol Error Rate');
    xlabel('SNR (dB)');
    ylabel('SER (%)');
    
    figure
    plot(snr_cnt, rmsevm_bpsk);
    hold on 
    plot(snr_cnt, rmsevm_qpsk);
    hold on 
    plot(snr_cnt, rmsevm_16qam);
    hold on 
    plot(snr_cnt, rmsevm_64qam);
    legend('BPSK','QPSK','16QAM','64QAM');
    title('Error Vector Magnitude');
    xlabel('SNR (dB)');
    ylabel('EVM');

    figure
    plot(snr_cnt, pkt_loss_bpsk);
    hold on
    plot(snr_cnt, pkt_loss_qpsk);
    hold on
    plot(snr_cnt, pkt_loss_16qam);
    hold on
    plot(snr_cnt, pkt_loss_64qam);
    legend('BPSK','QPSK','16QAM','64QAM');
    title('Packet Loss');
    xlabel('SNR (dB)');
    ylabel('Packet Loss (%)');
    

    toc
    

end


function [tx_vec, tx_syms]=transmit(tx_data, mod_scheme)
	TX_POWER_SCALE=1;
    switch mod_scheme
		case 'BPSK'
			modvec_bpsk = [-1 1];
			mod_fcn_bpsk = @(x) complex(modvec_bpsk(1+x),0);
			tx_syms = arrayfun(mod_fcn_bpsk, tx_data);
		case 'QPSK'
			modvec_qpsk = (1/sqrt(2)) .* [-1 1];
			mod_fcn_qpsk = @(x) complex(modvec_qpsk(1+bitshift(x, -1)), modvec_qpsk(1+mod(x, 2)));
			tx_syms = arrayfun(mod_fcn_qpsk, tx_data);
		case '16QAM'
			modvec_16qam = (1/sqrt(10)) .* [-3 -1 +3 +1];
			mod_fcn_16qam = @(x) complex(modvec_16qam(1+bitshift(x, -2)), modvec_16qam(1+mod(x,4)));
			tx_syms = arrayfun(mod_fcn_16qam, tx_data);
		case '64QAM'
			modvec_64qam = (1/sqrt(43)) .* [-7 -5 -1 -3 +7 +5 +1 +3];
			mod_fcn_64qam = @(x) complex(modvec_64qam(1+bitshift(x, -3)), modvec_64qam(1+mod(x,8)));
			tx_syms = arrayfun(mod_fcn_64qam, tx_data);
	end
	
	tx_vec=ifft(tx_syms);
	tx_vec=TX_POWER_SCALE .* tx_vec ./ max(abs(tx_vec));

end


function [rx_data, rx_syms]=receive(rx_vec,tx_syms,mod_scheme)
    TX_POWER_SCALE=1;
    rx_syms = fft(rx_vec);
    rs_syms= TX_POWER_SCALE .* rx_syms ./ max(abs(rx_syms));
    rx_syms = rx_syms ./ mean(abs(rx_syms)) .* mean(abs(tx_syms));
	switch mod_scheme
		case 'BPSK'
			demod_fcn_bpsk = @(x) double(real(x)>0);
			rx_data= arrayfun(demod_fcn_bpsk, rx_syms);
		case 'QPSK'
			demod_fcn_qpsk = @(x) double(2*(real(x)>0) + 1*(imag(x)>0));
			rx_data = arrayfun(demod_fcn_qpsk, rx_syms);    
		case '16QAM'
			demod_fcn_16qam = @(x) (8*(real(x)>0)) + (4*(abs(real(x))<0.6325)) + (2*(imag(x)>0)) + (1*(abs(imag(x))<0.6325));
			rx_data = arrayfun(demod_fcn_16qam, rx_syms);
		case '64QAM'
			demod_fcn_64qam = @(x) (32*(real(x)>0)) + (16*(abs(real(x))<0.6172)) + (8*((abs(real(x))<(0.9258))&&((abs(real(x))>(0.3086))))) + (4*(imag(x)>0)) + (2*(abs(imag(x))<0.6172)) + (1*((abs(imag(x))<(0.9258))&&((abs(imag(x))>(0.3086)))));
			rx_data = arrayfun(demod_fcn_64qam, rx_syms);
    end
end		

function [tx_data]=data_gen(mod_scheme, pkt_size)
    	bits_per_symbol=0;
	% Generating TX Data 
	switch mod_scheme
		case 'BPSK'
			bits_per_symbol=1;
		case 'QPSK'
			bits_per_symbol=2;
		case '16QAM'
			bits_per_symbol=4;
		case '64QAM'
			bits_per_symbol=6;
	end
	
	num_symbol=(pkt_size*8)/bits_per_symbol;
	tx_data = randi(2^bits_per_symbol, 1, num_symbol) - 1;
end

function [rx_vec, tx_vec_plot, rx_vec_plot]=channel_noise(tx_vec, snr)
    rx_vec = awgn(tx_vec,snr);	

    tx_vec_plot=tx_vec;
    x_val=1:30000;
    tx_vec_plot(numel(x_val)) = 0;
    rx_vec_plot=awgn(tx_vec_plot,snr);
end

function [sym_err_rate, bit_err_rate]= error_calc(tx_data, rx_data)
    [m,n]=size(rx_data);
    sym_err_rate = (sum(tx_data ~= rx_data)/(m*n))*100;

	rx_data_in_bit=decimalToBinaryVector(rx_data);
	tx_data_in_bit=decimalToBinaryVector(tx_data);
	bitErrors=0;
	[p,q]=size(rx_data_in_bit);
		for i=1:p
			for j=1:q
				if(rx_data_in_bit(i,j)~=tx_data_in_bit(i,j))
				bitErrors=bitErrors+1;
				end
			end
		end
	bit_err_rate = (bitErrors/(p*q))*100;
end

function plot_main(tx_vec, rx_vec, tx_syms, rx_syms)
    
    tx_vec_plot=tx_vec;
    x_val=1:30000;
    tx_vec_plot(numel(x_val)) = 0;

    % Plotting transmit waveform
    I_tx = real(tx_vec_plot); 
    Q_tx = imag(tx_vec_plot);
    figure;
    subplot(2,1,1);
    plot( gca, I_tx, 'b');
    title('TX - Inphase Component');
    
    subplot(2,1,2);
    plot( gca, Q_tx , 'r');
    title('TX - Quadrature Component');
    
    figure;
    subplot(2,1,1);
    plot( gca, I_tx(1:200), 'b');
    title('TX - Zoomed-in Inphase Component');
    
    subplot(2,1,2);
    plot( gca, Q_tx(1:200), 'r');
    title('TX - Zoomed-in Quadrature Component');


    % Plotting Received waveform
    I_rx = real(rx_vec); 
    Q_rx = imag(rx_vec);
    figure;
    subplot(2,1,1);
    plot( gca, I_rx, 'b');
    title('RX - Inphase Component');
    
    subplot(2,1,2);
    plot( gca, Q_rx , 'r');
    title('RX - Quadrature Component');
    
    figure;
    subplot(2,1,1);
    plot( gca, I_rx(1:200), 'b');
    title('RX - Zoomed-in Inphase Component');
    
    subplot(2,1,2);
    plot( gca, Q_rx(1:200), 'r');
    title('RX - Zoomed-in Quadrature Component');


    %Constellation diagram
    sPlotFig = scatterplot(rx_syms,1,0,'g+');
    hold on
    scatterplot(tx_syms,1,0,'rO',sPlotFig)
end

function [tx_data, tx_vec, tx_syms, rx_vec, rx_data, rx_syms, sym_err_rate, bit_err_rate, pkt_loss, tx_vec_plot, rx_vec_plot]=tx_and_rx(num_pkt, mod_scheme, pkt_size, snr)
    
    pkt_loss=0;
    
    for i=1:num_pkt

        % Generating Data
        [tx_data]=data_gen(mod_scheme, pkt_size);
	    
        % Applying modulation techniques
	    [tx_vec, tx_syms]=transmit(tx_data, mod_scheme);
    
	    % Adding noise
        [rx_vec, tx_vec_plot, rx_vec_plot]=channel_noise(tx_vec, snr);
	    
        % Demodulating and receiving the data at the receiver.
	    [rx_data, rx_syms] = receive(rx_vec,tx_syms, mod_scheme);
        
        % Calculating Errors
        [sym_err_rate_temp(i), bit_err_rate_temp(i)]= error_calc(tx_data, rx_data);
        
        % Calculating packet loss
        if(sym_err_rate_temp(i)~=0)
            pkt_loss=pkt_loss+1;
        end

    end
    sym_err_rate=mean(sym_err_rate_temp);
    bit_err_rate=mean(bit_err_rate_temp);
end
	