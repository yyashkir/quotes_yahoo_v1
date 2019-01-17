# 'Copyright © 2018 YASHKIR CONSULTING'
from pandas_datareader import data
import time
import datetime
import tkinter
from tkinter import *
from tkinter import messagebox
import matplotlib.pyplot as pyplot
from matplotlib.dates import date2num
#from matplotlib.finance import date2num
import numpy
import math
import os
import glob
import os.path
import shutil
import webbrowser
import warnings

def calc_date_span(date_from,date_to,t):
    t0 = datetime.datetime.strptime(date_from, '%Y-%m-%d')
    t1 = datetime.datetime.strptime(date_to, '%Y-%m-%d')
    for i in range(1,len(t)):
        if t[i-1]<=t0 and t[i]>t0:
            start_point = i-1
        if  t[i-1] ==t1 and t[i]>t1:
            end_point = i-1
        if  t[i-1] < t1 and t[i]>= t1:
            end_point = i
    return(start_point,end_point)


def save_2D_list(header,alist,file):
    f = open(file,'w')
    if len(header) != 0:
        f.write(header[0])
        for k in range(1, len(header)):
            f.write(',' + str(header[k]))
    for j in range(len(alist)):
        f.write('\n')
        f.write(alist[j][0])
        for k in range(1,len(alist[j])):
            f.write(','+ str(alist[j][k]))
    f.close()


def get_historical_prices(tickers,source,start_date,end_date):
    failed = ''
    done = ''
    new_data = []

    for i in range(len(tickers)):  # for every ticker
        res = None
        for s in range(3):
            try:
                p = data.DataReader(tickers[i], source,start_date,end_date)
                prices = list(p['Close'])
                dates = p.index.strftime('%Y-%m-%d').tolist()
                if i == 0:
                    new_data.append(dates)
                    new_data.append(prices)
                else:
                    iline = []                      # array of prices
                    # print('\n',tickers[i],'newdata ',new_data[0])
                    for ii in range(len(new_data[0])):  # for each date in dates for tickers[0]
                        # print(ii,'dates ',dates)
                        if new_data[0][ii] in dates:    # we check if this date is in dates for ticker[i]
                            ip = dates.index(new_data[0][ii])   # if yes then append corresponding price
                            iline.append(prices[ip])
                            # print('ip iline ',ip,iline)
                        else:
                            iline.append('NA')                  # or append NA
                            # print(ip,'ip NA iline ', iline)
                    new_data.append(iline)
                res = True
                if s == 0:
                    att_n = str(s+1)+' st'
                if s == 1:
                    att_n = str(s+1)+' nd'
                if s == 2:
                    att_n = str(s+1)+' rd'
                done = done + tickers[i] + '\ton '+att_n+' attempt\n'
                break
            except:
                pass

        if res == None:
            print('no price for:  ' + tickers[i],s)
            failed = failed + tickers[i] + '\n'
            if i == 0:
                messagebox.showinfo('ERROR', 'Data download failed for\n' + tickers[i] + '\nRestart the program')
                quit()
            u = len(dates)
            prices_NA = ['NA'] * u
            new_data.append(prices_NA)
    # for ii in range(len(tickers)):
    #     print(new_data[ii])
    if failed == '':
        messagebox.showinfo('Download summary', 'Prices downloaded for:\n'+done)
    else:
        messagebox.showinfo('Download summary',  'Prices downloaded for:\n'+done + '\n'+ 'Failed download for: '+failed)
    # new_data[0] has all dates, new_data[1,2,...] have all prices (one line for a ticker)
    # transposing list of lists
    V = []
    for j in range(len(new_data[0])):
        U = []
        for i in range(len(new_data)):
            try:
                U.append(new_data[i][j])
            except:
                U.append('NA')
        V.append(U)
    new_data = V
    return(new_data)

def get_header_tickers_names(stock_list):
    tickers = []  # as in menu
    names = []  # as in menu
    header = ['Date']
    for k in range(len(stock_list)):
        tickers.append(stock_list[k][0])  #
        names.append(stock_list[k][1])
        header.append(tickers[k] + '_' + names[k])
    return (header,tickers, names)

def  datafile_update(stock_list_in_file,quotes_in_file, stock_list_menu,new_data,quote_file_name):
   
    # stock_list_in_file:   [['NYSE:TD', 'Toronto Dominion Bank'], ['GDP', 'Good
    # quotes_in_file[0]:    ['2016-11-19', '47.26', '5.49', '54.2', '769.2']
    # stock_list_menu:      [['NYSE:TD', 'Toronto Dominion Bank'], ['NYSE:BNS', 'Bank of Nova Scotia'], ['GOOG', 'Google']]
    # new_data[0]):         ['2016-11-25', 47.689999999999998, 53.939999999999998, 761.67999999999995
    tics_rem = []
    tics_add = []
    numb_of_stocks_in_file = len(stock_list_in_file)
    # print('numb stocks in file before rem/add: ',numb_of_stocks_in_file)
    for k in range(len(stock_list_in_file)):
        if stock_list_in_file[k] not in stock_list_menu:
            tics_rem.append(stock_list_in_file[k])
    for k in range(len(stock_list_menu)):
        if stock_list_menu[k] not in stock_list_in_file:
            tics_add.append(stock_list_menu[k])
    # print('remove ',tics_rem,'  add ',tics_add)
    for p in range(len(tics_rem) - 1, -1, -1):      # removing columns (tickers and prices)
        k = stock_list_in_file.index(tics_rem[p])
        # print('removing ',k,stock_list_in_file[k])
        for l in quotes_in_file:
            l.pop(k + 1)  # k+1 because element 0 is a date (removing price in a row "l" in column "k+1"
        stock_list_in_file.pop(k)  # shortened (stock in column "k" removed
        # print('numb stocks after rem/add: ', len(stock_list_in_file))
        # print('numb columns in quotes table: ',len(quotes_in_file[0]))
    quotes_updated = []
    for j in range(len(quotes_in_file)):
        row = []
        for k in range(len(stock_list_in_file) + 1):  # number of columns is now different
            row.append(quotes_in_file[j][k])  # building a row of prices
        quotes_updated.append(row)  # quotes_updated building up
    stock_list_updated = stock_list_in_file[
                         :]  # stock_list_in_file is actually updated list (not as in original file), renaming...
 
    for k in range(len(tics_add)):
        stock_list_updated.append(tics_add[k])

    for j in range(len(quotes_updated)):    # addidng columns for tics_add data
        for s in range(len(tics_add)):
            quotes_updated[j].append('NA')

    # no new data downloaded, tickers with corresponding columns removed and/or added (NA)
    if len(new_data)==0:
        header, tickers, names = get_header_tickers_names(stock_list_updated)
        save_2D_list(header, quotes_updated, quote_file_name)
        return (stock_list_updated, quotes_updated, tickers, names)

    # adding lines with additional dates and no data ('NA')
    start_date_in_file = datetime.datetime.strptime(quotes_updated[0][0], "%Y-%m-%d").date()
    end_date_in_file = datetime.datetime.strptime(quotes_updated[-1][0], "%Y-%m-%d").date()
    start_date_in_new_data = new_data[0][0]
    start_block = []
    for f in range(len(new_data)):
        date_in_new_data = datetime.datetime.strptime(new_data[f][0], "%Y-%m-%d").date()
        if date_in_new_data > end_date_in_file:
            new_line = []
            new_line.append(new_data[f][0])
            for k in range(1,len(quotes_updated[0])):
                new_line.append('NA')
            quotes_updated.append(new_line)
        if date_in_new_data < start_date_in_file:
            new_line = []
            new_line.append(new_data[f][0])
            for k in range(1,len(quotes_updated[0])):
                new_line.append('NA')
            start_block.append(new_line)
    _all = []
    for j in range(len(start_block)):
        _all.append(start_block[j])
    for j in range(len(quotes_updated)):
        _all.append(quotes_updated[j])
    quotes_updated = _all

    # adding data to quotes_updated from new_data
    dates_new = []
    for f in range(len(new_data)):
        dates_new.append(new_data[f][0])

    for j in range(len(quotes_updated)):
        tj = quotes_updated[j][0]
        if tj in dates_new:
            s = dates_new.index(tj)
            for k in range(1, len(stock_list_menu)+1):
                if quotes_updated[j][k] == 'NA' or quotes_updated[j][k]=='nan' and new_data[s][k] != 'NA':
                    quotes_updated[j][k] = new_data[s][k]
    header, tickers, names = get_header_tickers_names(stock_list_menu)
    save_2D_list(header, quotes_updated, quote_file_name)
    return(stock_list_menu,quotes_updated,tickers,names)

def select(lstbox):
    reslist = list()
    seleccion = lstbox.curselection()
    for i in seleccion:
        entrada = lstbox.get(i)
        reslist.append(entrada)
    return(reslist)

def show_last_prices(fname):
    f = open(fname, "r")
    txt = f.read()
    split_txt = txt.split('\n')
    f.close()
    clean_txt = []
    for j in range(len(split_txt)):
        if len(split_txt[j]) != 0:
            clean_txt.append(split_txt[j])
    header = clean_txt[0].split(',')
    prices = clean_txt[-1].split(',')
    txt = str(prices[0])+'\nPrice\tEquity'
    for k in range(1,len(header)):
        try:
            txt = txt + '\n' + str("{0:.2f}".format(float(prices[k]))) + '\t' + header[k]
        except:
            txt = txt + '\n' + str(prices[k]) + '\t' + header[k]

    root = Tk()
    root.title('The last date prices')
    txt_win_height = len(header) + 2
    content = Text(root, width=60, height=txt_win_height, wrap='none', font=("Times", 10), bg="light cyan", fg="black")
    content.insert(INSERT, txt)
    content.grid(row=0, column=0, sticky=N + S)
    yscroll = Scrollbar(root, command=content.yview, orient=VERTICAL)
    yscroll.grid(row=0, column=1, rowspan=txt_win_height, sticky=N + S)
    content.configure(yscrollcommand=yscroll.set)

    xscroll = Scrollbar(root, command=content.xview, orient=HORIZONTAL)
    xscroll.grid(row=1, column=0, columnspan=100, sticky=W + E)
    content.configure(xscrollcommand=xscroll.set)
    Button(root, bg='black', fg='red', text='Close window ',
           font=("Arial", 12, "bold"), command=lambda: root.destroy()).grid(row=3, column=0, columnspan=2, sticky=W + E)
    root.mainloop()

def trade_on_MA(tics_to_display,tickers,names,delt,quote_list,date_from,date_to,buy_criterium,sell_criterium,mov_aver_days_short,mov_aver_days_long,trend_days):
    year_days = 365
    week   = ['Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday']
    nm = len(quote_list)
    date_now = datetime.datetime.strptime(quote_list[nm-1][0], "%Y-%m-%d").date()
    wdn = date_now.weekday()
    weekday = week[wdn]
    down = float(buy_criterium) / 100.
    up = float(sell_criterium) / 100.

    for j in range(len(tickers)):      # for each ticker ...
        if tickers[j] not in tics_to_display:
            continue
        moving_short = []
        moving_long = []
        y = []                         # array of prices
        t = []
        t_float = []
        v0 = 0
        for i in range(nm):           # along time line...
            if type(quote_list[i][j+1]) != float:
                print(tickers[j] + ': no price on ' + quote_list[i][0])
                continue                            # to skip points where price is not available

            y.append(quote_list[i][j+1])     # from data list quote_list
            t.append(datetime.datetime.strptime(quote_list[i][0], '%Y-%m-%d'))
            t_float.append(date2num(datetime.datetime.strptime(quote_list[i][0], '%Y-%m-%d')))

            mov_aver_days_short_modif = min(len(y),mov_aver_days_short)
            mov_aver_days_long_modif = min(len(y),mov_aver_days_long)

            mov_aver_short = sum(y[len(y)-mov_aver_days_short_modif: len(y)]) / mov_aver_days_short_modif
            mov_aver_long = sum(y[len(y)-mov_aver_days_long_modif: len(y)]) / mov_aver_days_long_modif

            moving_short.append(mov_aver_short)
            moving_long.append(mov_aver_long)

        # checking dates
        if datetime.datetime.strptime(date_from, '%Y-%m-%d') < t[0]:
            date_from = t[0].strftime('%Y-%m-%d')
        if datetime.datetime.strptime(date_to, '%Y-%m-%d') > t[-1]:
            date_to = t[-1].strftime('%Y-%m-%d')
        if  datetime.datetime.strptime(date_to, '%Y-%m-%d') < t[0]      or \
            datetime.datetime.strptime(date_from, '%Y-%m-%d') > t[-1]   or \
            datetime.datetime.strptime(date_from, '%Y-%m-%d') >= datetime.datetime.strptime(date_to, '%Y-%m-%d'):
            message = 'The date range for ' + tickers[j] + ' is wrong. \nNo graph, sorry...\nTo continue?'

            if messagebox.askquestion("Continue/Break", message) == 'yes':
                continue
            else:
                break

        start_point,end_point = calc_date_span(date_from,date_to,t)

        t_float = t_float[start_point:end_point + 1]
        t0 = float(t_float[0])
        for s in range(len(t_float)):
            t_float[s] = float(t_float[s]) - t0    #relative time in days starting at zero

        t = t[start_point:end_point+1]
        y = y[start_point:end_point+1]
        moving_short = moving_short[start_point:end_point+1]
        moving_long = moving_long[start_point:end_point+1]
        start_point = 0
        end_point = len(t) - 1

        trade_value = []
        buy_schedule = []
        sell_schedule = []
        trade_schedule = []
        buy_points = []
        sell_points = []
        cash = y[0]
        stock_num = 0
        trade_value.append(y[0])
        trade_schedule.append(t[0])
        trade_log = '\n ' + tickers[j]
        for i in range(trend_days,len(t)):
            # in thi s loop the end point is "i"
            # trend calc starts here
            start_point_trend = max(0,i - trend_days)
            t_date_trend = t[start_point_trend:i + 1]
            t_trend = t_float[start_point_trend:i + 1]
            y_trend = y[start_point_trend:i + 1]

            y_lin = []
            try:
                z = numpy.polyfit(t_trend, y_trend, 1)
                p = numpy.poly1d(z)
                for ii in range(len(t_trend)):
                    y_lin.append(p(t_trend[ii]))
            except:
                print('no lin fit for ' + tickers[j])
            slope = p.coeffs[0]

            #buy event
            if stock_num == 0 and slope >0 and y[i] < moving_short[i] * (1. - down):
                stock_num = cash / y[i]
                trade_schedule.append(t[i])
                trade_value.append(cash)
                buy_schedule.append(t[i].date())
                buy_points.append(y[i])
                trade_log = trade_log + '\n ' + str("{0:.2f}".format(stock_num)) + ' stocks bought at ' + str("{0:.2f}".format(y[i])) + ' on ' + str(t[i].date())

            #sell event
            if stock_num >0 and stock_num * y[i] >= trade_value[-1] * (1. + up):
                trade_log = trade_log + '\n ' + str("{0:.2f}".format(stock_num)) + ' stocks sold at ' + str("{0:.2f}".format(y[i])) + ' on ' + str(t[i].date())
                cash = stock_num * y[i]
                trade_schedule.append(t[i])
                trade_value.append(cash)
                stock_num = 0
                sell_schedule.append(t[i].date())
                sell_points.append(y[i])

        # print(trade_log)

        trade_schedule.append(t[-1])
        if stock_num == 0:
            trade_value.append(cash)
        else:
            trade_value.append(stock_num * y[-1])

        days_span = t_float[-1] - t_float[0]
        year_span = days_span / year_days
        gain = trade_value[-1] / trade_value[0] - 1

        if gain > -1:
            gain_annualized = "{0:.2f}".format(math.log(1 + gain) / year_span * 100)
        else:
            gain_annualized = ''
        gain_pc = "{0:.2f}".format(( trade_value[-1] / trade_value[0] - 1 )* 100)

        graphtext = 'Lines: Price (grey) / MA '+str(mov_aver_days_short)+' d (red)\n' \
                    'Circles: bought (blue) / sold (red); Portfolio (black dots)'
        trade_info = str(len(trade_schedule))+' trades \nGain '+str(gain_pc)+' % (' + str(gain_annualized) + ' % - annualized)'

        fig1 = pyplot.figure(1)
        period_y= int(days_span / year_days)
        period_extra_days = int(days_span - period_y * year_days)
        period_txt = str(period_y)+' years + ' + str(period_extra_days) + ' days'
        date_info = str(t[0].date()) + ' to ' + str(t[-1].date()) + ' (' + period_txt+')'
        pyplot.title('Backtesting: '+tickers[j]+' ('+names[j]+') \n'+ date_info, color='navy',fontweight='bold')
        pyplot.ylabel('Price',color='navy',fontweight='bold',rotation = 45)
        pyplot.xlabel(graphtext, color='navy', fontweight='bold', rotation=0)
        ax = pyplot.gca()
        x0 = 0.05
        y0 = 0.95
        hor = 'left'
        pyplot.text(x0,y0,trade_info,horizontalalignment= hor,verticalalignment='top',transform = ax.transAxes,color='blue') #,fontweight='')
        pyplot.grid()
        graph = fig1.add_subplot(111)
        graph.patch.set_facecolor('white')
        graph.tick_params(axis='x', colors='blue')
        graph.tick_params(axis='y', colors='black')

        graph.plot_date(t,y,'-',markersize=4,color='grey')
        graph.plot_date(t,moving_short,'-',color='red')
        # graph.plot_date(t,moving_long,'-',color='black')
        graph.plot_date(buy_schedule, buy_points, 'o', color='blue', fillstyle='none')
        graph.plot_date(sell_schedule, sell_points, 'o', color = 'red', fillstyle='none')
        graph.plot_date(trade_schedule, trade_value, 'o',color = 'black', fillstyle = 'full')
        pyplot.gcf().autofmt_xdate()

        fig_hight = len(trade_schedule) * 0.25
        fig2 = pyplot.figure(2,figsize=(6,fig_hight))
        graph2 = fig2.add_subplot(111)
        graph2.plot([0,2],[0,1],alpha = 0.0)
        pyplot.xticks([0,2]," ")
        pyplot.yticks([0,1], " ")
        pyplot.text(0, 1, trade_log, horizontalalignment='left', verticalalignment='top',
                    color='blue')
        pyplot.title('Backtesting: ' + tickers[j] + ' (' + names[j] + ') \n' + date_info, color='navy',
                     fontweight='bold')
        pyplot.show()

        img_filename = tickers[j].replace(':','-') + '.png'
        pyplot.savefig(img_filename)
        result = messagebox.askokcancel(title=None,message='Next chart ?')
        if result == FALSE :
            pyplot.close('all')
            break

def trade_linear(tics_to_display,tickers,names,delt,quote_list,date_from,date_to,buy_criterium,sell_criterium,mov_aver_days_short,mov_aver_days_long,trend_days):
    year_days = 365
    week   = ['Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday']
    nm = len(quote_list)
    date_now = datetime.datetime.strptime(quote_list[nm-1][0], "%Y-%m-%d").date()
    wdn = date_now.weekday()
    weekday = week[wdn]
    down = float(buy_criterium) / 100.
    up = float(sell_criterium) / 100.

    for j in range(len(tickers)):      # for each ticker ...
        if tickers[j] not in tics_to_display:
            continue
        moving_short = []
        moving_long = []
        y = []                         # array of prices
        t = []
        t_float = []
        v0 = 0
        for i in range(nm):           # along time line...
            if type(quote_list[i][j+1]) != float:
                print(tickers[j] + ': no price on ' + quote_list[i][0])
                continue                            # to skip points where price is not available

            y.append(quote_list[i][j+1])     # from data list quote_list
            t.append(datetime.datetime.strptime(quote_list[i][0], '%Y-%m-%d'))
            t_float.append(date2num(datetime.datetime.strptime(quote_list[i][0], '%Y-%m-%d')))

            mov_aver_days_short_modif = min(len(y),mov_aver_days_short)
            mov_aver_days_long_modif = min(len(y),mov_aver_days_long)

            mov_aver_short = sum(y[len(y)-mov_aver_days_short_modif: len(y)]) / mov_aver_days_short_modif
            mov_aver_long = sum(y[len(y)-mov_aver_days_long_modif: len(y)]) / mov_aver_days_long_modif

            moving_short.append(mov_aver_short)
            moving_long.append(mov_aver_long)

        # checking dates
        if datetime.datetime.strptime(date_from, '%Y-%m-%d') < t[0]:
            date_from = t[0].strftime('%Y-%m-%d')
        if datetime.datetime.strptime(date_to, '%Y-%m-%d') > t[-1]:
            date_to = t[-1].strftime('%Y-%m-%d')
        if  datetime.datetime.strptime(date_to, '%Y-%m-%d') < t[0]      or \
            datetime.datetime.strptime(date_from, '%Y-%m-%d') > t[-1]   or \
            datetime.datetime.strptime(date_from, '%Y-%m-%d') >= datetime.datetime.strptime(date_to, '%Y-%m-%d'):
            message = 'The date range for ' + tickers[j] + ' is wrong. \nNo graph, sorry...\nTo continue?'
            if messagebox.askquestion("Continue/Break", message) == 'yes':
                continue
            else:
                break

        start_point,end_point = calc_date_span(date_from,date_to,t)

        t_float = t_float[start_point:end_point + 1]
        t0 = float(t_float[0])
        for s in range(len(t_float)):
            t_float[s] = float(t_float[s]) - t0    #relative time in days starting at zero

        t = t[start_point:end_point+1]
        y = y[start_point:end_point+1]
        moving_short = moving_short[start_point:end_point+1]
        moving_long = moving_long[start_point:end_point+1]
        start_point = 0
        end_point = len(t) - 1

        trade_value = []
        buy_schedule = []
        sell_schedule = []
        trade_schedule = []
        buy_points = []
        sell_points = []
        cash = y[0]
        stock_num = 0
        trade_value.append(y[0])
        trade_schedule.append(t[0])
        trade_log = '\n ' + tickers[j]
        for i in range(trend_days,len(t)):
            # in thi s loop the end point is "i"
            # trend calc starts here
            start_point_trend = max(0,i - trend_days)
            t_date_trend = t[start_point_trend:i + 1]
            t_trend = t_float[start_point_trend:i + 1]
            y_trend = y[start_point_trend:i + 1]

            y_lin = []
            try:
                z = numpy.polyfit(t_trend, y_trend, 1)
                p = numpy.poly1d(z)
                for ii in range(len(t_trend)):
                    y_lin.append(p(t_trend[ii]))
            except:
                print('no lin fit for ' + tickers[j])
            slope = p.coeffs[0]

            upper_lin = y_lin[-1] * (1. + up)
            lower_lin = y_lin[-1] * (1. - down)

            #buy event
            if stock_num == 0 and slope >0 and y[i] < lower_lin:
                stock_num = cash / y[i]
                trade_log = trade_log + '\n' + str("{0:.2f}".format(stock_num)) + ' stocks bought at ' + str("{0:.2f}".format(y[i])) + ' on ' + str(t[i].date())
                trade_schedule.append(t[i])
                trade_value.append(cash)
                buy_schedule.append(t[i].date())
                buy_points.append(y[i])

            #sell event
            if stock_num > 0 and y[i] >= upper_lin and slope > 0:
                trade_log = trade_log + '\n' + str("{0:.2f}".format(stock_num)) + ' stocks sold at ' + str("{0:.2f}".format(y[i])) + ' on ' + str(t[i].date())
                cash = stock_num * y[i]
                trade_schedule.append(t[i])
                trade_value.append(cash)
                stock_num = 0
                sell_schedule.append(t[i].date())
                sell_points.append(y[i])
        print(trade_log)
        trade_schedule.append(t[-1])
        if stock_num == 0:
            trade_value.append(cash)
        else:
            trade_value.append(stock_num * y[-1])



        days_span = t_float[-1] - t_float[0]
        year_span = days_span / year_days
        gain = trade_value[-1] / trade_value[0] - 1

        if gain > -1:
            gain_annualized = "{0:.2f}".format(math.log(1 + gain) / year_span * 100)
        else:
            gain_annualized = ''
        gain_pc = "{0:.2f}".format(( trade_value[-1] / trade_value[0] - 1 )* 100)

        graphtext = 'Lines: Price (grey) / Moving Average '+str(mov_aver_days_short)+' days (red) \n' \
                    'Circles: bought (blue) / sold (red);  Portfolio (black)'
        trade_info = str(len(trade_schedule))+' trades \n Gain '+str(gain_pc)+' % (' + str(gain_annualized) + ' % - annualized)'

        fig1 = pyplot.figure(1)
        period_y= int(days_span / year_days)
        period_extra_days = int(days_span - period_y * year_days)
        period_txt = str(period_y)+' years + ' + str(period_extra_days) + ' days'
        date_info = str(t[0].date()) + ' to ' + str(t[-1].date()) + ' (' + period_txt+')'
        pyplot.title('Backtesting: '+tickers[j]+' ('+names[j]+') \n'+ date_info, color='navy',fontweight='bold')
        pyplot.ylabel('Price',color='navy',fontweight='bold',rotation = 45)
        pyplot.xlabel(graphtext, color='navy', fontweight='bold', rotation=0)
        ax = pyplot.gca()
        x0 = 0.05
        y0 = 0.95
        hor = 'left'
        pyplot.text(x0,y0,trade_info,horizontalalignment= hor,verticalalignment='top',transform = ax.transAxes,color='blue') #,fontweight='')
        pyplot.grid()
        graph = fig1.add_subplot(111)
        graph.patch.set_facecolor('white')
        graph.tick_params(axis='x', colors='blue')
        graph.tick_params(axis='y', colors='black')

        graph.plot_date(t,y,'-',markersize=4,color='grey')
        graph.plot_date(t,moving_short,'-',color='red')
        # graph.plot_date(t,moving_long,'-',color='black')
        graph.plot_date(buy_schedule, buy_points, 'o', color='blue', fillstyle='none')
        graph.plot_date(sell_schedule, sell_points, 'o', color = 'red', fillstyle='none')
        graph.plot_date(trade_schedule, trade_value, 'o',color = 'black', fillstyle = 'full')
        pyplot.gcf().autofmt_xdate()
        fig_hight = len(trade_schedule) * 0.3
        fig2 = pyplot.figure(2,figsize=(6,fig_hight))
        graph2 = fig2.add_subplot(111)
        graph2.plot([0,2],[0,1],alpha = 0.0)
        pyplot.xticks([0,2]," ")
        pyplot.yticks([0,1], " ")
        pyplot.text(0, 1, trade_log, horizontalalignment='left', verticalalignment='top',
                    color='blue')
        pyplot.title('Backtesting: ' + tickers[j] + ' (' + names[j] + ') \n' + date_info, color='navy',
                     fontweight='bold')
        pyplot.show()

        img_filename = tickers[j].replace(':','-') + '.png'
        pyplot.savefig(img_filename)
        result = messagebox.askokcancel(title=None,message='Next chart ?')
        if result == FALSE :
            pyplot.close('all')
            break

def trade_updown(tics_to_display,tickers,names,delt,quote_list,date_from,date_to,buy_criterium,sell_criterium,mov_aver_days_short,mov_aver_days_long,trend_days):
    year_days = 365
    week   = ['Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday']
    nm = len(quote_list)
    date_now = datetime.datetime.strptime(quote_list[nm-1][0], "%Y-%m-%d").date()
    wdn = date_now.weekday()
    weekday = week[wdn]
    down = float(buy_criterium) / 100.
    up = float(sell_criterium) / 100.

    for j in range(len(tickers)):      # for each ticker ...
        if tickers[j] not in tics_to_display:
            continue
        moving_short = []
        moving_long = []
        y = []                         # array of prices
        t = []
        t_float = []
        v0 = 0
        for i in range(nm):           # along time line...
            if type(quote_list[i][j+1]) != float:
                print(tickers[j] + ': no price on ' + quote_list[i][0])
                continue                            # to skip points where price is not available

            y.append(quote_list[i][j+1])     # from data list quote_list
            t.append(datetime.datetime.strptime(quote_list[i][0], '%Y-%m-%d'))
            t_float.append(date2num(datetime.datetime.strptime(quote_list[i][0], '%Y-%m-%d')))

            mov_aver_days_short_modif = min(len(y),mov_aver_days_short)
            mov_aver_days_long_modif = min(len(y),mov_aver_days_long)

            mov_aver_short = sum(y[len(y)-mov_aver_days_short_modif: len(y)]) / mov_aver_days_short_modif
            mov_aver_long = sum(y[len(y)-mov_aver_days_long_modif: len(y)]) / mov_aver_days_long_modif

            moving_short.append(mov_aver_short)
            moving_long.append(mov_aver_long)

        # checking dates
        if datetime.datetime.strptime(date_from, '%Y-%m-%d') < t[0]:
            date_from = t[0].strftime('%Y-%m-%d')
        if datetime.datetime.strptime(date_to, '%Y-%m-%d') > t[-1]:
            date_to = t[-1].strftime('%Y-%m-%d')
        if  datetime.datetime.strptime(date_to, '%Y-%m-%d') < t[0]      or \
            datetime.datetime.strptime(date_from, '%Y-%m-%d') > t[-1]   or \
            datetime.datetime.strptime(date_from, '%Y-%m-%d') >= datetime.datetime.strptime(date_to, '%Y-%m-%d'):
            message = 'The date range for ' + tickers[j] + ' is wrong. \nNo graph, sorry...\nTo continue?'
            if messagebox.askquestion("Continue/Break", message) == 'yes':
                continue
            else:
                break

        start_point,end_point = calc_date_span(date_from,date_to,t)

        t_float = t_float[start_point:end_point + 1]
        t0 = float(t_float[0])
        for s in range(len(t_float)):
            t_float[s] = float(t_float[s]) - t0    #relative time in days starting at zero

        t = t[start_point:end_point+1]
        y = y[start_point:end_point+1]
        moving_short = moving_short[start_point:end_point+1]
        moving_long = moving_long[start_point:end_point+1]
        start_point = 0
        end_point = len(t) - 1

        trade_value = []
        buy_schedule = []
        sell_schedule = []
        trade_schedule = []
        buy_points = []
        sell_points = []
        cash = y[0]
        stock_num = 1
        trade_value.append(y[0])
        trade_schedule.append(t[0])
        trade_price = y[0]
        trade_log =  '\n ' + tickers[j]
        trade_log = trade_log + '\n ' + str("{0:.2f}".format(stock_num)) + ' stocks bought at ' + str("{0:.2f}".format(y[0])) + ' on ' + str(t[0].date())
        for i in range(trend_days,len(t)):
            # in thi s loop the end point is "i"
            # trend calc starts here
            start_point_trend = max(0,i - trend_days)
            t_date_trend = t[start_point_trend:i + 1]
            t_trend = t_float[start_point_trend:i + 1]
            y_trend = y[start_point_trend:i + 1]

            y_lin = []
            try:
                z = numpy.polyfit(t_trend, y_trend, 1)
                p = numpy.poly1d(z)
                for ii in range(len(t_trend)):
                    y_lin.append(p(t_trend[ii]))
            except:
                print('no lin fit for ' + tickers[j])
            slope = p.coeffs[0]

            upper_lin = y_lin[-1] * (1. + up)
            lower_lin = y_lin[-1] * (1. - down)

            # print(i,y[i],trade_price)
            #buy event
            if stock_num == 0 and y[i] < trade_price * (1. - down):
                trade_price = y[i]
                stock_num = cash / y[i]
                trade_schedule.append(t[i])
                trade_value.append(cash)
                buy_schedule.append(t[i].date())
                buy_points.append(y[i])
                trade_log = trade_log + '\n ' + str("{0:.2f}".format(stock_num)) + ' stocks bought at ' + str("{0:.2f}".format(y[i])) + ' on ' + str(t[i].date())
                # print(str(t[i].date()) + ': Bought '+str("{0:.2f}".format(stock_num))+' stocks at '+str("{0:.2f}".format(y[i])))

            #sell event
            if stock_num > 0 and y[i] >= trade_price * (1. + up) and slope <=0:
                trade_log = trade_log + '\n ' + str("{0:.2f}".format(stock_num)) + ' stocks sold at ' + str("{0:.2f}".format(y[i])) + ' on ' + str(t[i].date())
                # print(str(t[i].date()) + ': Sold ' + str("{0:.2f}".format(stock_num)) + ' stocks at ' + str(
                #     "{0:.2f}".format(y[i])))
                trade_price = y[i]
                cash = stock_num * y[i]
                trade_schedule.append(t[i])
                trade_value.append(cash)
                stock_num = 0
                sell_schedule.append(t[i].date())
                sell_points.append(y[i])
        print(trade_log)
        print('\n')

        trade_schedule.append(t[-1])
        if stock_num == 0:
            trade_value.append(cash)
        else:
            trade_value.append(stock_num * y[-1])

        days_span = t_float[-1] - t_float[0]
        year_span = days_span / year_days
        gain = trade_value[-1] / trade_value[0] - 1

        if gain > -1:
            gain_annualized = "{0:.2f}".format(math.log(1 + gain) / year_span * 100)
        else:
            gain_annualized = ''
        gain_pc = "{0:.2f}".format(( trade_value[-1] / trade_value[0] - 1 )* 100)

        graphtext = 'Lines: Price (grey) / Moving Average '+str(mov_aver_days_short)+' days (red)\n' \
                    'Circles: bought (blue) / sold (red); Portfolio (black)'
        trade_info = str(len(trade_schedule))+' trades \n Gain '+str(gain_pc)+' % (' + str(gain_annualized) + ' % - annualized)'
        fig1 = pyplot.figure(1)
        period_y= int(days_span / year_days)
        period_extra_days = int(days_span - period_y * year_days)
        period_txt = str(period_y)+' years + ' + str(period_extra_days) + ' days'
        date_info = str(t[0].date()) + ' to ' + str(t[-1].date()) + ' (' + period_txt+')'
        pyplot.title('Backtesting: '+tickers[j]+' ('+names[j]+') \n'+ date_info, color='navy',fontweight='bold')
        pyplot.ylabel('Price',color='navy',fontweight='bold',rotation = 45)
        pyplot.xlabel(graphtext, color='navy', fontweight='bold', rotation=0)
        ax = pyplot.gca()
        x0 = 0.05
        y0 = 0.95
        hor = 'left'
        pyplot.text(x0,y0,trade_info,horizontalalignment= hor,verticalalignment='top',transform = ax.transAxes,color='blue') #,fontweight='')
        pyplot.grid()
        graph = fig1.add_subplot(111)
        graph.patch.set_facecolor('white')
        graph.tick_params(axis='x', colors='blue')
        graph.tick_params(axis='y', colors='black')

        graph.plot_date(t,y,'-',markersize=4,color='grey')
        graph.plot_date(t,moving_short,'-',color='red')
        # graph.plot_date(t,moving_long,'-',color='black')
        graph.plot_date(buy_schedule, buy_points, 'o', color='blue', fillstyle='none')
        graph.plot_date(sell_schedule, sell_points, 'o', color = 'red', fillstyle='none')
        graph.plot_date(trade_schedule, trade_value, 'o',color = 'black', fillstyle = 'full')
        pyplot.gcf().autofmt_xdate()

        fig_hight = len(trade_schedule) * 0.25
        fig2 = pyplot.figure(2,figsize=(6,fig_hight))
        graph2 = fig2.add_subplot(111)
        graph2.plot([0,2],[0,1],alpha = 0.0)
        pyplot.xticks([0,2]," ")
        pyplot.yticks([0,1], " ")
        pyplot.text(0, 1, trade_log, horizontalalignment='left', verticalalignment='top',
                    color='blue')
        pyplot.title('Backtesting: ' + tickers[j] + ' (' + names[j] + ') \n' + date_info, color='navy',
                     fontweight='bold')
        pyplot.show()

        img_filename = tickers[j].replace(':','-') + '.png'
        pyplot.savefig(img_filename)
        result = messagebox.askokcancel(title=None,message='Next chart ?')
        if result == FALSE :
            pyplot.close('all')
            break

def monitor(tics_to_display,tickers,names,delt,quote_list,date_from,date_to,buy_criterium,sell_criterium,mov_aver_days_short,mov_aver_days_long,trend_days):
    year_days = 365
    week   = ['Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday']
    nm = len(quote_list)
    date_now = datetime.datetime.strptime(quote_list[nm-1][0], "%Y-%m-%d").date()
    wdn = date_now.weekday()
    weekday = week[wdn]
    down = float(buy_criterium) / 100.
    up = float(sell_criterium) / 100.
    date1 = date_from
    date2 = date_to
    for j in range(len(tickers)):      # for each ticker ...
        date_from = date1   #restoring date span for next ticker
        date_to = date2
        if tickers[j] not in tics_to_display:
            continue
        moving_short = []
        moving_long = []
        y = []                         # array of prices
        t = []
        t_float = []
        v0 = 0
        for i in range(nm):           # along time line...
            if type(quote_list[i][j+1]) != float:
                print(tickers[j] + ': no price on ' + quote_list[i][0])
                continue                            # to skip points where price is not available

            y.append(quote_list[i][j+1])     # from data list quote_list
            t.append(datetime.datetime.strptime(quote_list[i][0], '%Y-%m-%d'))
            t_float.append(date2num(datetime.datetime.strptime(quote_list[i][0], '%Y-%m-%d')))

            mov_aver_days_short_modif = min(len(y),mov_aver_days_short)
            mov_aver_days_long_modif = min(len(y),mov_aver_days_long)

            mov_aver_short = sum(y[len(y)-mov_aver_days_short_modif: len(y)]) / mov_aver_days_short_modif
            mov_aver_long = sum(y[len(y)-mov_aver_days_long_modif: len(y)]) / mov_aver_days_long_modif

            moving_short.append(mov_aver_short)
            moving_long.append(mov_aver_long)
        if datetime.datetime.strptime(date_from, '%Y-%m-%d') < t[0]:
            date_from = t[0].strftime('%Y-%m-%d')
        if datetime.datetime.strptime(date_to, '%Y-%m-%d') > t[-1]:
            date_to = t[-1].strftime('%Y-%m-%d')
        if  datetime.datetime.strptime(date_to, '%Y-%m-%d') < t[0]      or \
            datetime.datetime.strptime(date_from, '%Y-%m-%d') > t[-1]   or \
            datetime.datetime.strptime(date_from, '%Y-%m-%d') >= datetime.datetime.strptime(date_to, '%Y-%m-%d'):
            message = 'The date range for ' + tickers[j] + ' is wrong. \nNo graph, sorry...\nTo continue?'
            if messagebox.askquestion("Continue/Break", message) == 'yes':
                continue
            else:
                break

        start_point,end_point = calc_date_span(date_from,date_to,t)

        t_float = t_float[start_point:end_point + 1]
        t0 = float(t_float[0])
        for s in range(len(t_float)):
            t_float[s] = float(t_float[s]) - t0    #relative time in days starting at zero

        t = t[start_point:end_point+1]
        y = y[start_point:end_point+1]
        moving_short = moving_short[start_point:end_point+1]
        moving_long = moving_long[start_point:end_point+1]
        start_point = 0
        end_point = len(t) - 1

        start_point_trend = max(0,end_point - trend_days)
        t_date_trend = t[start_point_trend:end_point + 1]
        t_trend = t_float[start_point_trend:end_point + 1]
        y_trend = y[start_point_trend:end_point + 1]

        y_lin = []
        try:
            z = numpy.polyfit(t_trend, y_trend, 1)
            p = numpy.poly1d(z)
            for ii in range(len(t_trend)):
                y_lin.append(p(t_trend[ii]))
        except:
            print('no lin fit for ' + tickers[j])
        slope = p.coeffs[0]

# to buy or to sell ?
        # by Moving Average criterium
        hint = 'MA: Wait'
        #buy event
        if slope >0 and y[end_point] < moving_short[end_point] * (1. - down):
            hint = 'MA: Buy'
        #sell event
        if slope <= 0 and y[end_point] > moving_short[end_point] * (1. + up):
            hint = 'MA: Sell'
        buy_sell = hint

        # by Linear approximation criterium
        hint = '\nLinear: Wait'
        if slope > 0 and y[-1] < y_lin[-1] * (1 - down):
            hint = '\nLinear: Buy'
        if slope > 0 and y[-1] >= y_lin[-1] * (1 + up):
            hint = '\nLinear: Sell'
        buy_sell = buy_sell + hint

        # by up/down criterium
        hint = '\nUp/Down: Wait'
        if y[-1] < y[0] * (1. - down):
            hint = '\nUp/Down: Buy'
        if y[-1] > y[0] * (1. + up):
            hint = '\nUp/Down: Sell'
        buy_sell = buy_sell + hint

        days_span = t_float[-1] - t_float[0]
        year_span = days_span / year_days

        last_price = str("{0:.2f}".format(float(y[-1])))
        graphtext = 'Lines: Moving Average: '+str(mov_aver_days_short)+' d (red) / '+str(mov_aver_days_long)+' d (blue), Linear regression (dashed)\n'+ 'Circles: Price (grey)'
        trade_info = buy_sell
        fig = pyplot.figure()
        period_y= int(days_span / year_days)
        period_extra_days = int(days_span - period_y * year_days)
        period_txt = str(period_y)+' years + ' + str(period_extra_days) + ' days'
        date_info = str(t[0].date()) + ' to ' + str(t[-1].date()) + ' (' + period_txt+')'
        pyplot.title(tickers[j]+' ('+names[j]+'): $'+ last_price +  '\n'+ date_info, color='navy',fontweight='bold')
        pyplot.ylabel('Price',color='navy',fontweight='bold',rotation = 45)
        pyplot.xlabel(graphtext, color='navy', fontweight='bold', rotation=0)
        ax = pyplot.gca()
        x0 = 0.05
        y0 = 0.95
        hor = 'left'
        pyplot.text(x0,y0,trade_info,horizontalalignment= hor,verticalalignment='top',transform = ax.transAxes,color='blue') #,fontweight='')
        pyplot.grid()
        graph = fig.add_subplot(111)
        graph.patch.set_facecolor('white')
        graph.tick_params(axis='x', colors='blue')
        graph.tick_params(axis='y', colors='black')

        graph.plot_date(t,y,'o',markersize=4,color='grey')
        graph.plot_date(t,moving_short,'-',color='red')
        graph.plot_date(t,moving_long,'-',color='blue')
        graph.plot_date(t_date_trend, y_lin, '--', color='black')

        pyplot.gcf().autofmt_xdate()
        fig.show()
        img_filename = tickers[j].replace(':','-') + '.png'
        pyplot.savefig(img_filename)
        result = messagebox.askokcancel(title=None,message='Next chart ?')
        if result == FALSE :
            pyplot.close('all')
            break

def run(lstbox,valores,pars,tickers_titles,quotes_in_file,quote_file_name,stock_list_in_file):
# stock_list_in_file =    [['NYSE:TD', 'Toronto Dominion Bank'], ['GDP', 'Goodrich Petroleum Corp'],...]
    task = pars[0].get()
    date_from = pars[1].get()
    date_to = pars[2].get()
    buy_criterium = pars[3].get()
    sell_criterium = pars[4].get()
    mov_aver_days_short = int(pars[5].get())      #
    mov_aver_days_long = int(pars[6].get())      #               #
    backup_file_name = str(pars[7].get())
    mode = str(pars[8].get())
    trend_days = int(pars[9].get())

    if quote_file_name != backup_file_name:
        shutil.copyfile(quote_file_name, backup_file_name)
    if quote_file_name == backup_file_name:
        if messagebox.askyesno('Warning','To overwrite '+'quotes.csv ?'):
            shutil.copyfile(quote_file_name, backup_file_name)
        else:
            pass
    stocks_text_menu = tickers_titles.get("1.0", END) # stock list from menu (could be different from stock (from file)
                    # NYSE:TD_Toronto Dominion Bank
                    # GDP_Goodrich Petroleum Corp
                    # NYSE:BNS_Bank of Nova Scotia
                    # GOOG_Google
    try:
        start_date = datetime.datetime.strptime(date_from, '%Y-%m-%d')
    except:
        messagebox.showinfo('Date Error','The start date '+str(date_from)+' format is wrong')
        return(0)
    try:
        end_date = datetime.datetime.strptime(date_to, '%Y-%m-%d')
    except:
        messagebox.showinfo('Date Error','The end date '+str(date_to)+' format is wrong')
        return(0)

    date = time.strftime("%d/%m/%Y")   # getting present date
    tics_names = stocks_text_menu.split('\n')

    stock_list_menu = []
    tickers = []
    for k in range(len(tics_names)):
        tic_nam = tics_names[k].split('_')
        if len(tic_nam)>1:
            stock_list_menu.append([tic_nam[0],tic_nam[1]])
            tickers.append(tic_nam[0])

    new_data = []
    if task == 'Yes':    # for stocks in menu
        source = 'yahoo'    #'google'
        new_data = get_historical_prices(tickers, source,start_date,end_date)

    stock_list_updated, quotes_updated,tickers,names = datafile_update(stock_list_in_file,quotes_in_file, stock_list_menu,new_data,quote_file_name)
    tics_2choos_from = tickers[0]
    for k in range(1,len(tickers)):
        tics_2choos_from = tics_2choos_from + ' '+tickers[k]
    valores.set(tics_2choos_from)
    tics_to_display = select(lstbox)
    dates = []
    delt = []
    date0 = datetime.datetime.strptime(quotes_updated[0][0], "%Y-%m-%d").date() # the oldest date in file
    for i in range(0,len(quotes_updated)):
        if(not quotes_updated[i][0]==''):
            dates.append(quotes_updated[i][0])          # loading all dates to 'dates'
            date_now = datetime.datetime.strptime(quotes_updated[i][0], "%Y-%m-%d").date()   # converting data string to data format
            delt.append((date_now - date0).days)                                   # difference in days (converting to days)
            for j in range(1,len(quotes_updated[i])):
                try:
                    quotes_updated[i][j] = float(quotes_updated[i][j])      # converting string to float
                except ValueError:
                    # print( quotes_updated[i][j])
                    quotes_updated[i][j] = 'NA'

    if mode == 'Monitor':
        monitor(tics_to_display, tickers, names, delt, quotes_updated, date_from, date_to, buy_criterium,
                    sell_criterium, mov_aver_days_short, mov_aver_days_long, trend_days)
        # make_graphs(tics_to_display,tickers,names,delt,quotes_updated,date_from,date_to,buy_criterium,sell_criterium,mov_aver_days_short,mov_aver_days_long,trend_days)        # graphs

    if mode == 'Trade LINEAR':
        trade_linear(tics_to_display, tickers, names, delt, quotes_updated, date_from, date_to, buy_criterium, sell_criterium,
                 mov_aver_days_short, mov_aver_days_long, trend_days)

    if mode == 'Trade UPDOWN':
        trade_updown(tics_to_display, tickers, names, delt, quotes_updated, date_from, date_to, buy_criterium,
                     sell_criterium,
                     mov_aver_days_short, mov_aver_days_long, trend_days)

    if mode == 'Trade MA':
        trade_on_MA(tics_to_display,tickers,names,delt,quotes_updated,date_from,date_to,buy_criterium,sell_criterium,mov_aver_days_short,mov_aver_days_long,trend_days)
    messagebox.showinfo("showinfo",'\nYou may continue or quit')
    return(0)

def read_menu_conf(mfile):
    try:
        f = open(mfile, "r")  # opening
    except IOError:
        tkinter.messagebox.showinfo("showinfo", 'The file ' + mfile + ' is missing')
        quit()
    names = []
    values = []
    while 'TRUE':
        line = f.readline()
        if not line:
            break
        tokens = line.strip().split(':')
        names.append(tokens[0])
        values.append(tokens[1])
    values_list = []
    for k in range(len(values)):
        values_list.append(values[k].split(','))
    f.close()
    return(names,values_list)

def read_quotes(qfile):
    try:
        f = open(qfile, "r")  # opening historical data
    except IOError:
        tkinter.messagebox.showinfo("showinfo", 'The file ' + qfile + ' is missing')
        quit()
    line0 = f.readline()
    tic_name = line0.strip().split(',')
    tic_name.pop(0)
    stock_list_in_file = []
    for k in range(len(tic_name)):
        stock_list_in_file.append(tic_name[k].split('_'))
    quotes_in_file = []
    while 'TRUE':
        line = f.readline()
        if not line:
            break
        if (len(line) > 1):
            tokens = line.strip().split(',')  # csv are split to list
            quotes_in_file.append(tokens)  # loading historicalprices to quotes_in_file
    f.close()
    return(stock_list_in_file,quotes_in_file)

def save_widget_text(content,fname):
    text = content.get("1.0", END)
    f = open(fname, "w")
    f.write(text)
    f.close()

def show_edit_txt_file(file_name,title,width,image_file):
    f = open(file_name, "r")
    txt = f.read()
    f.close()
    root = Tk()
    root.title(title)
    txt_win_height = 60
    content = Text(root, width=width, height=txt_win_height, wrap='none', font=("Times", 10), bg="light cyan", fg="black")
    content.insert(INSERT, txt)
    content.grid(row=0, column=0, sticky=N + S)
    yscroll = Scrollbar(root, command=content.yview, orient=VERTICAL)
    yscroll.grid(row=0, column=1, rowspan=1, sticky=N + S)
    content.configure(yscrollcommand=yscroll.set)
    xscroll = Scrollbar(root, command=content.xview, orient=HORIZONTAL)
    xscroll.grid(row=1, column=0, columnspan=1, sticky=W + E)
    content.configure(xscrollcommand=xscroll.set)
    if image_file != '':
        photo = PhotoImage(master=root, file='yy.gif')
        content.image_create(END, image=photo, align='baseline')
    Button(root, bg='black', fg='navy', text='Save', font=("Arial", 12, "bold"), command=lambda: save_widget_text(content, file_name)).grid(row=3, column=0, columnspan=2, sticky=E)
    Button(root, bg='black', fg='red', text='Close window (if saved the program restart is required)',font=("Arial", 12, "bold"), command=lambda: root.destroy()).grid(row=4, column=0, columnspan=2, sticky=W + E)
    root.mainloop()

def exit(man):
    for f in glob.glob("*.png"):
        if os.path.isfile(f):
            os.remove(f)
            print(f+' removed')
    man.destroy()

def callback():
    webbrowser.open_new(r"http://www.yashkir.com")

def menu_function(file_range,menu_names,menu_values, stock_list_in_file,quotes_in_file,quote_file_name):
    file_info = 'Data file: '+ quote_file_name + '. Date range: from '+ file_range[0] + '  to  '+ file_range[1]
    man = Tk()
    pars = []
    k = 0
    for j in range(len(menu_names)):
        if (len(menu_values[j])== 1):
            Label(man,text=menu_names[j],borderwidth=1, font=("Arial", 10, "bold")).grid(row=k,column=0,sticky=E)
            pars.append(StringVar())
            pars[k].set(menu_values[j][0])
            Entry(man,textvariable= pars[k],bd = 2,justify=LEFT).grid(row=k,column=1,sticky=E+W)
            k = k + 1
        if (len(menu_values[j]) > 1):  # menu items with a set of values in the drop-down list
            Label(man, text=menu_names[j], borderwidth=1, font=("Arial", 10, "bold")).grid(row=k, column=0,sticky=E)
            pars.append(StringVar())
            types = menu_values[j]
            pars[k].set(types[0])
            w = OptionMenu(man,pars[k], *types)
            w.config(width=15,bg='light blue',activebackground='red')
            w.grid(row=k, column=1)
            k = k + 1

    stock_text = ''
    tickers = []
    for j in range(len(stock_list_in_file)):   # [['NYSE:TD', 'Toronto Dominion Bank'], ['GDP', 'Goodrich Petr ....]]
        if j==0:
            stock_text = stock_text + str(stock_list_in_file[j][0]) + '_' + str(stock_list_in_file[j][1])
        else:
            stock_text = stock_text + '\n' + str(stock_list_in_file[j][0]) + '_' + str(stock_list_in_file[j][1])
        tickers.append(stock_list_in_file[j][0])
    k = k + 1

    font_b = ("Arial", 10, "bold")
    man.title('Market quotes')
    txt_win_height = 30 # len(stock_list_in_file)  #max(50,len(stock_list_in_file))
    tickers_titles = Text(man, width=40, height=txt_win_height, font=("Times", 12, 'italic'),background="light cyan", foreground="blue")
    tickers_titles.insert(INSERT, stock_text)
    tickers_titles.grid(row=1, column=2, rowspan=20,sticky= N+S)
    Label(man, text="Ticker_Company Name\n(Remove or add if desired)", borderwidth=1, background="white",foreground="navy")\
        .grid(row=0, column=2, sticky=N+W + E)

    tics_in_file = ''
    for k in range(len(stock_list_in_file)):
        tics_in_file = tics_in_file + ' ' + stock_list_in_file[k][0]

    # selection menu starts here
    valores = StringVar()
    valores.set(tics_in_file)
    lstbox = Listbox(man, listvariable=valores, selectmode=MULTIPLE,  height=30,background="yellow", foreground="blue",font=("Arial", 10, "bold"))
    for k in range(len(tics_in_file)):
        lstbox.select_set(k)
    lstbox.grid(column=3, row=1,rowspan=20,sticky= N+S)

    Label(man, text="Unselect/Select tickers \nfor displaying charts", borderwidth=1,background="yellow", foreground="navy")\
        .grid(row=0, column=3, sticky=W+E)

    Button(man,height=1,bg='light blue',activebackground='red', fg='dark blue',bd = 2,text='RUN',font=("Arial",12,"bold"),command = lambda:
        run(lstbox,valores,pars,tickers_titles,quotes_in_file,quote_file_name,stock_list_in_file))\
        .grid(row=k+2,column=0,columnspan=2,sticky= W+E)
    # row k+3   is not used
    Button(man, height=1, bg='blue',activebackground='red', fg='navy', text='HELP', font=("Arial", 10, "bold"),command=lambda:
    show_edit_txt_file('help', 'Help',100,'yy.gif')).grid(row=k+4, column=3, sticky=W + E)
    Label(man, text=file_info, borderwidth=5, bg="steel blue",activebackground='red', fg="navy")\
        .grid(row=k + 2, column=2, columnspan=1,sticky=W + E)
    Button(man, height=1, bg='gray30',activebackground='red', fg='navy', text='Open/edit data file',
           font=("Arial", 10, "bold"), command=lambda: show_edit_txt_file(quote_file_name, 'Data', 200,'')).grid(row=k + 4, column=0, columnspan=1, sticky=W + E)
    Button(man, height=1, bg='black',activebackground='red', fg='navy', text='Last date prices',
           font=("Arial", 10, "bold"), command=lambda: show_last_prices(quote_file_name))\
        .grid(row=k + 4, column=1,columnspan=1, sticky=W + E)
    Button(man,bg='blue',activebackground='red',fg='red',text='Quit',font=("Arial",10,"bold"),command = lambda:
        exit(man)).grid(row=k+2,column=3,columnspan=1,sticky=W+ E)
    Button(man, height=1, bg='gray30',activebackground='red', fg='navy', text='Copyright © 2018 YASHKIR CONSULTING', font=("Arial", 10, "bold"),command=lambda:
        callback()).grid(row=k+4, column=2, columnspan = 1, sticky=W+E)

    man.mainloop()

# starting!
menu_file = 'menu.conf'
quote_file_name = 'quotes.csv'
menu_names,menu_values = read_menu_conf(menu_file)
stock_list_in_file,quotes_in_file = read_quotes(quote_file_name)
file_range = [quotes_in_file[0][0],quotes_in_file[-1][0]]
ind = menu_names.index('From')
menu_values[ind][0] = quotes_in_file[-1][0]
ind = menu_names.index('To')
now = datetime.date.today().strftime('%Y-%m-%d')
menu_values[ind][0] = now

menu_function(file_range,menu_names,menu_values, stock_list_in_file,quotes_in_file,quote_file_name)

# 'Copyright © 2017 YASHKIR CONSULTING'
#   Last updated 25/01/2018
